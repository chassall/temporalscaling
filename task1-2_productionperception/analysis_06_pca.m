% Run PCA on residual EEG
% Also does figures - may move this to a separate script
%
% Project: Temporal Scaling
% Other m-files required: EEGLAB, subtightplot.m, makefigure.m

% Author: Cameron Hassall, Department of psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% March 2020; Last revision: 12-Jul-2022

close all; clear all; clc; rng(41); % for reproducibility

% Set data and results folders (change as needed)
dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task1-2_productionperception/data';
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

% Analysis settings: 'full', 'fixed', or 'scaled'
whichAnalysis = 'scaled'; % 'full' EEG, 'fixed' only (scaled regressed out), 'scaled' only (fixed regressed out)

participants = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};
nParticipants = length(participants);
whichTask = 'production';
whichMethod = 'box';
srate = 200;
behavLimits = [0.2 5];
aveChannel = 0;
targetTimes = [0.8 1.65 2.5];
trialLengths = [0 0.8; 0 1.65; 0 2.5];

% Axis settings
axisSettings.fontSize = 8;
axisSettings.fontName = 'Arial';
margin = 0;
mmeanwindow = 25;

%% Quantile settings

% Channels: P3, CP5, CP1 (24, 18, 19)

% pcaChannels = [9 10 11 15 19 20 21]; % Channels to submit to the PCA
% pcaChannels = [24, 18, 19];
pcaChannelStrings = {'Fp1'	'Fp2'	'F7'	'F3'	'Fz'	'F4'	'F8'	'FC5'	'FC1'	'FCz'	'FC2'	'FC6'	'T7'	'C3'	'Cz'	'C4'	'T8'	'CP5'	'CP1'	'CPz'	'CP2'	'CP6'	'P7'	'P3'	'Pz'	'P4'	'P8'	'POz'	'O1'	'Oz'	'O2'	'AFz'};

if aveChannel
    whichChannel = 1;
else
    whichChannelString = 'P3'; % Channel of interest
    whichChannel = find(strcmp(pcaChannelStrings,whichChannelString));
end


% if aveChannel
%     whichChannel = 1;
% else
%     % whichChannel = 4; % Channel of interest 10: fcz, 15: Cz
%     whichChannel = 1; 
% end
whichPC = 2; % PC of interest
nQuantiles = 2; % Number of quantile boundaries
quantileNames = {'q1','q2','q3'};
withinNames = [{'1'};{'2'};{'3'}];
rmModel = 'q1-q3 ~ 1';
numPCs = 6;

conditionStrings = {'Short','Medium','Long'};
triggers = {{12},{22},{32}};
for c = 1:3
    trialTimes{c} = trialLengths(c,1):1/srate:trialLengths(c,2)-1/srate;
end
baseline = [1000*margin*(1/srate) 80];

allDataForPCA{1} = [];
allDataForPCA{2} = [];
allDataForPCA{3} = [];

allPCAScores = [];
allPCACoeff{1} = [];
allPCACoeff{2} = [];
allPCACoeff{3} = [];

allPCAScores2{1} = [];
allPCAScores2{2} = [];
allPCAScores2{3} = [];

allBinAverages = {};

explainedByParticipant = [];

allRs = [];
allMeanRTs = [];
allSRTs = [];
allMRTs = [];
allLRTs = [];
for iParticipant = 1:nParticipants
    disp(['participant ' participants{iParticipant}]);
    subjString = ['sub-' participants{iParticipant}];

    % Load preprocessed data
    preprocessedFolder = fullfile(dataFolder,'derivatives/eegprep/',subjString);
    % preprocessedFolder = ['./data/derivatives/eegprep/sub-' participants{iParticipant}];
    preprocessedFile = [subjString '_task-temporalscaling_eegprep.mat'];
    disp(['loading ' preprocessedFile]);
    load(fullfile(preprocessedFolder,preprocessedFile), 'EEG');
    % Remove ocular components using the results of ICLabel
    eyeLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Eye'));
    [~,I] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
    whichOnes = find(I == eyeLabel); % Max out of all possible labels
    EEG = pop_subcomp(EEG,whichOnes,0);
    
    % Load the residuals, which will be in EEGLAB format (called rEEG)
    residualsFolder = fullfile(dataFolder,'derivatives/eegresidual/',subjString);
    % residualsFolder = ['./data/derivatives/eegresidual/sub-' participants{iParticipant}];
    residualsFile = [subjString '_task-' whichTask '_eegresidual' whichMethod];
    load(fullfile(residualsFolder,residualsFile),'rEEGWithTS','rEEGWithFixed','toExclude');
    
    switch whichAnalysis
        case 'full'
            rEEG = EEG;
        case 'fixed'
            rEEG = rEEGWithFixed;
        case 'scaled'
            rEEG = rEEGWithTS;
    end
    
    pcaChannels = eeg_chaninds(rEEG,pcaChannelStrings);

    % Remove bad trials
    badShortTrials = [];
    sCount = 1;
    badMediumTrials = [];
    mCount = 1;
    badLongTrials = [];
    lCount = 1;
    for i = 1:length(rEEG.event)
        if ismember(rEEG.event(i).type,[12 22 32])
            startTime = rEEG.event(i).latency;
            endTime = rEEG.event(i+1).latency;
            thesePoints = startTime:endTime;
            switch rEEG.event(i).type
                case 12
                    if any(ismember(thesePoints,find(toExclude)))
                        badShortTrials = [badShortTrials sCount];
                    end
                    sCount = sCount + 1;
                case 22
                    if any(ismember(thesePoints,find(toExclude)))
                        badMediumTrials = [badMediumTrials mCount];
                    end
                    mCount = mCount + 1;
                case 32
                    if any(ismember(thesePoints,find(toExclude)))
                        badLongTrials = [badLongTrials lCount];
                    end
                    lCount = lCount + 1;
            end
        end
    end
    
    % Tried a filter 19 April 2020
    % rEEG = pop_eegfiltnew(rEEG, [], filterCutoff);
    
    % View residuals?
    % pop_select(rEEG);
    
    % Load the behavioural data
    behFolder = fullfile(dataFolder,subjString,'beh');
    % behFolder = ['./data/rawdata/sub-' participants{iParticipant} '/beh'];
    behFile = [subjString '_task-temporalscaling_beh.tsv'];
    % thisData = load(fullfile(behFolder,behFile));
    opts = delimitedTextImportOptions("NumVariables", 14);
    opts.Delimiter = "\t";
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    thisData = readtable(fullfile(behFolder,behFile), opts);
    thisData = table2array(thisData); % Do this because we originally worked with an array

    
    % Get RT for each condition (short, medium, long)
    behArtifacts = [];
    rts = [];
    for c = 1:3
        rts(:,c) = thisData(thisData(:,1) == c, 4);
        
        for i = 1:length(thisData(thisData(:,1) == c, 4))
            
            thisRT = rts(i,c);
            % if thisRT < targetTimes(c)-behavMargin(c) || thisRT > targetTimes(c)+behavMargin(c)
            if thisRT < behavLimits(1) || thisRT > behavLimits(2)
                behArtifacts(i,c) = 1;
            else
                behArtifacts(i,c) = 0;
            end
        end
    end
    
    
    % Grab residual EEG for each condition and trial
    % Also detect artifacts
    eEEG = {};
    eegArtifacts = [];
    for c = 1:3
        eEEG{c} = pop_epoch(rEEG,triggers{c},trialLengths(c,:));
        eEEG{c} = pop_rmbase(eEEG{c},baseline);
        eegArtifacts(:,c) = find_artifacts(eEEG{c});
    end
    
    % Define any artifact as either behavioural or EEG
    isAnyArtifact = behArtifacts | eegArtifacts;
    
    % Make bins a different way
    cleanSlowRTs = rts(~isAnyArtifact(:,1),1);
    cleanMediumRTs = rts(~isAnyArtifact(:,2),2);
    cleanFastRTs = rts(~isAnyArtifact(:,3),3);
    
    % New as of 2021-02-10
    allMeanRTs(iParticipant,:) = [mean(cleanSlowRTs) mean(cleanMediumRTs) mean(cleanFastRTs)];
    allSRTs = [allSRTs cleanSlowRTs'];
    allMRTs = [allMRTs cleanMediumRTs'];
    allLRTs = [allLRTs cleanFastRTs'];
    
    if nQuantiles == 1
        slowQuantiles = quantile(cleanSlowRTs,0.5);
        mediumQuantiles = quantile(cleanMediumRTs,0.5);
        fastQuantiles = quantile(cleanFastRTs,0.5);
    else
        slowQuantiles = quantile(cleanSlowRTs,nQuantiles);
        mediumQuantiles = quantile(cleanMediumRTs,nQuantiles);
        fastQuantiles = quantile(cleanFastRTs,nQuantiles);
    end
    
    % * * * Do categories a different way * * *
    slowCategories = nan(size(cleanSlowRTs));
    mediumCategories = nan(size(cleanMediumRTs));
    fastCategories = nan(size(cleanFastRTs));
    for r = 1:length(cleanSlowRTs)
        slowCategories(r) = sum(find(cleanSlowRTs(r) <= [slowQuantiles max(cleanSlowRTs)],1));
    end
    for r = 1:length(cleanMediumRTs)
        mediumCategories(r) = sum(find(cleanMediumRTs(r) <= [mediumQuantiles max(cleanMediumRTs)],1));
    end
    for r = 1:length(cleanFastRTs)
        fastCategories(r) = sum(find(cleanFastRTs(r) <= [fastQuantiles max(cleanFastRTs)],1));
    end
    
    % theseCategories should be zero for behavioural or EEG artifacts, and
    % a number from 1-nQuantiles+1 otherwise
    theseCategories = zeros(size(rts));
    theseCategories(~isAnyArtifact(:,1),1) = slowCategories;
    theseCategories(~isAnyArtifact(:,2),2) = mediumCategories;
    theseCategories(~isAnyArtifact(:,3),3) = fastCategories;
    
    % Make bin averages
    binAverages = {};
    binTotals = {};
    for b = 1:nQuantiles+1
        theseTrials = theseCategories == b;
        
        for c = 1:3
            thisEEG = eEEG{c};
            badTrials = find_artifacts(thisEEG);
            binAverages{c,b} = mean(thisEEG.data(pcaChannels,:,theseTrials(:,c)' & ~badTrials),3);
            allBinAverages{iParticipant,c,b} = mean(thisEEG.data(pcaChannels,:,theseTrials(:,c)' & ~badTrials),3);
            binCount{c,b} = sum(theseTrials(:,c));
        end
        
    end
    
    % Do PCA
    pcaScores = [];
    pcaCoeff = {};
    pcaScore = {};
    for c = 1:3
        
        pcaData = [];
        for b = 1:nQuantiles+1
            if aveChannel
                pcaData = [pcaData; mean(binAverages{c,b},1)];
            else
                pcaData = [pcaData; binAverages{c,b}];
            end
        end
        
        % For the PCA on all participant data
        allDataForPCA{c} = [allDataForPCA{c}; pcaData];
        
        % Do the PCA for this participant
        [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(pcaData,'Centered',false);
        explainedByParticipant(iParticipant,c,:) = EXPLAINED(1:9);
        
        % Flip the sign on the first PC as needed
        if mean(COEFF(:,1)) > 0
            COEFF(:,1) = -COEFF(:,1);
            SCORE(:,1) = -SCORE(:,1);
        end
        
        % Flip the sign on the second PC as needed
        if mean(COEFF(1:end/2,2)) > mean(COEFF(:,2))
            COEFF(:,2) = -COEFF(:,2);
            SCORE(:,2) = -SCORE(:,2);
        end
        
        % Store PCA scores for this condition
        if aveChannel
            pcaScores(c,:,:) =  reshape(SCORE(:,whichPC),1,[]);
        else
            pcaScores(c,:,:) =  reshape(SCORE(:,whichPC),length(pcaChannels),[]);
        end
        pcaCoeff{c} = COEFF;
        pcaScore{c} = SCORE;
    end
    allPCAScores(iParticipant,:,:,:) = pcaScores;
    for c = 1:3
        allPCACoeff{c}(iParticipant,:,:) = pcaCoeff{c};
        allPCAScores2{c}(iParticipant,:,:) = pcaScore{c};
    end
    
    % Do a PCA on all trials for this participant
    if 0
        figure();
        for c = 1:3
            tempChannel = 15;
            
            badTrialsBeh = behArtifacts(:,c);
            badTrialsEEG = find_artifacts(eEEG{c});
            thisPCAData = squeeze(eEEG{c}.data(tempChannel,:,:));
            thisPCAData(:,badTrialsEEG|badTrialsBeh) = [];
            thisPCAData = thisPCAData';
            [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(thisPCAData,'Centered',false);
            
            % Here I flip the sign on the second PC as needed
            %             if mean(COEFF(1:end/2,1)) < mean(COEFF(:,1))
            %                 COEFF(:,1) = -COEFF(:,1);
            %                 SCORE(:,1) = -SCORE(:,1);
            %             end
            %             if mean(COEFF(1:end/2,2)) > mean(COEFF(:,2))
            %                 COEFF(:,2) = -COEFF(:,2);
            %                 SCORE(:,2) = -SCORE(:,2);
            %             end
            if mean(COEFF(:,1)) > 0
                COEFF(:,1) = -COEFF(:,1);
                SCORE(:,1) = -SCORE(:,1);
            end
            if mean(COEFF(1:end/2,2)) > mean(COEFF(:,2))
                COEFF(:,2) = -COEFF(:,2);
                SCORE(:,2) = -SCORE(:,2);
            end
            
            theseTimes = rts;
            theseTimes(badTrialsEEG|badTrialsBeh) = [];
            theseBins = theseCategories(:,c);
            theseBins(badTrialsEEG|badTrialsBeh) = [];
            
            % Bin averages
            theseData = [];
            for b = 1:nQuantiles+1
                theseBinTimes(b) = mean(theseTimes(theseBins==b));
                theseBinScores(b) = mean(SCORE(theseBins==b,2));
                
                theseRecData = SCORE(theseBins==b,1) * COEFF(:,1)' + SCORE(theseBins==b,2) * COEFF(:,2)';
                
                theseBinCoeff(b) = mean(COEFF(theseBins==b,2));
                theseData(b,:) = mean(theseRecData,1);
            end
        end
    end
end

%% Mean Variance
tVal = abs(tinv(0.025,nParticipants-1));
meanVar = squeeze(mean(explainedByParticipant(:,:,1:2),1));
stdVar = squeeze(std(explainedByParticipant(:,:,1:2),[],1));
ciVar = tVal * stdVar ./ sqrt(nParticipants);
disp('Mean Variance Explained');
disp(meanVar);
disp(meanVar - ciVar);
disp(meanVar + ciVar);

%% ERPs
plot_settings;
for c = 1:3
    for b = 1:nQuantiles+1
        binEEG{c,b} = [];
    end
end

for iParticipant = 1:nParticipants
    
    for c = 1:3
        
        for b = 1:nQuantiles+1
            thisData = allBinAverages{iParticipant,c,b}(whichChannel,:);
            binEEG{c,b} = [binEEG{c,b}; thisData];
        end
    end
end

% ERPs
binERPs{1} = [];
binERPs{2} = [];
binERPs{3} = [];
for c = 1:3
    for b = 1:nQuantiles+1
        binERP{c,b} = mean(binEEG{c,b},1);
    end
end

subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.01 0.1], [0.3 0.1], [0.07 0.03]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis.
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins

plotColours = cbrewer('qual','Dark2',nQuantiles+1);
makefigure(19,8);
titles = {'Short','Medium','Long'};
for c = 1:3
    subplot(1,3,c);
    
    for q = 1:nQuantiles+1
        plot(trialTimes{c},binERP{c,q},'Color',plotColours(q,:),'LineWidth',plotLineWidth);
        hold on;
    end
    
    if c == 2
        l = legend('very early','early','on time','late','very late');
        l.Box = 'off';
        l.Orientation = 'horizontal';
        l.Position = [0.1978    0.0118    0.6418    0.0670];
    end
    
    ax = gca();
    ax.Box = 'off';
    % ax.Color = 'k';
    ax.FontSize = plotFontSize;
    ax.XAxis.Label.String = 'Time (ms)';
    ax.YAxis.Label.String = 'Voltage (\muv)';
    ax.Title.String = titles{c};
end

% print(['.tiffs/fig_' whichTask '_' whichAnalysis '_residual_erps_by_quantile.tiff'],'-dtiff','-r300');

%%
figure();
for c = 1:3
    for q = 1:nQuantiles+1
        plot(trialTimes{c},binERP{c,q},'LineWidth',plotLineWidth);
        hold on;
    end
end

%%
figure();
allData = [];
whichCondition = 3;
for iParticipant = 1:nParticipants
    thisData = allBinAverages{iParticipant,whichCondition,2}(whichChannel,:);
    allData = [allData; thisData];
end

%%
fig = makefigure(9,7);
sAve = [];
mAve = [];
lAve = [];
for q = 1:nQuantiles+1
    sAve = [sAve; binERP{1,q}];
    mAve = [mAve; binERP{2,q}];
    lAve = [lAve; binERP{3,q}];
end
sAve = mean(sAve,1);
mAve = mean(mAve,1);
lAve = mean(lAve,1);
plot_settings;
plot(trialTimes{1},sAve,'LineWidth',plotLineWidth,'Color',plotColours(1,:));
hold on;
plot(trialTimes{2},mAve,'LineWidth',plotLineWidth,'Color',plotColours(2,:));
plot(trialTimes{3},lAve,'LineWidth',plotLineWidth,'Color',plotColours(3,:));
ax = gca;
ax.Box = 'off';
ax.XAxis.Label.String = 'Time (s)';
ax.YAxis.Label.String = 'Voltage (\muV)';
l = legend('Short','Medium','Long');
l.Box = 'off';
l.Location = 'southwest';
% print(['./tiffs/fig_' whichTask '_' whichAnalysis '_residual_by_condition.tiff'],'-dtiff','-r300');

%% PCA on all participant data
allCoeff = {};
allScores = {};
allExplained = [];
for c = 1:3
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(allDataForPCA{c},'Centered',false);
    COEFF(:,1) = -COEFF(:,1);
    SCORE(:,1) = -SCORE(:,1);
    allExplained(:,c) = EXPLAINED(1:10);
    
    if mean(COEFF(1:end/2,2)) > mean(COEFF(:,2))
        COEFF(:,2) = -COEFF(:,2);
        SCORE(:,2) = -SCORE(:,2);
    end
    allCoeff{c} = COEFF;
    allScores{c} = SCORE;
end


%% Scree plot
plot_settings;
fig = makefigure(9,6);
for c = 1:3
    plot(allExplained(1:5,c),'LineWidth',plotLineWidth,'Color',plotColours(c,:));
    hold on;
end
l = legend('Short','Medium','Long');
l.Box = 'off';
ax = gca;
ax.Box = 'off';
ax.YLabel.String = 'Proportion explained (%)';
ax.XLabel.String = 'PC Number';

%% Visualize PCA (Manuscript Figure)
plot_settings
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.01 0.04], [0.14 0.1], [0.1 0.03]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis.
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins
%

pcaColours = cbrewer('div','BrBG',256);
pcaColours = flip(pcaColours);
makefigure(8.51,7.19);
plotIndices = {[1],[2 3],[4 5 6]};
for c = 1:3
    xs = trialTimes{c};
    ys = 1:size(allDataForPCA{c},1)/nParticipants;
    
    subplot(1,6,plotIndices{c});
    
    thisPColour = pcolor(xs,ys,allDataForPCA{c}(ys,:));
    thisPColour.LineStyle = 'none';
    % thisPColour.FaceColor = 'interp';
    caxis([-20 20]);
    ax = gca;
    set_axis(ax,axisSettings);
    ax.FontSize = plotFontSize;
    % ax.Colormap = pcaColours; % Stopped working after upgrading MATLAB :(
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'Quantiles/Sensors';
    ax.YTickLabel = {};
    if ismember(c,[2 3])
        ax.YAxis.Visible = 'off';
    end
end
cb = colorbar();
cb.Label.String = 'Voltage (\muV)';

%%
figure();
for c = 1:3
    disp([min(allScores{c}(:,2)) max(allScores{c}(:,2))]);
    
    
    subplot(1,3,c)
    plot(allCoeff{c}(:,1:3));
    legend();
end

%%
whichCoefficients = [];
for c = 1:3
    thisMin = min(allScores{c}(:,2));
    thisMax = max(allScores{c}(:,2));
    thisMin = -10;
    thisMax = 100;
    whichCoefficients = [whichCoefficients; linspace(thisMin,thisMax,5)];
end

figure();
plot_settings;
for c = 1:3
    subplot(1,3,c);
    plotColours = cbrewer('div','RdYlBu',size(whichCoefficients,2));
    plot(mean(allDataForPCA{c},1),'k:','LineWidth',plotLineWidth);
    hold on;
    for w = 1:size(whichCoefficients,2)
        theseScores = allScores{c};
        theseScores(:,2) = whichCoefficients(c,w);
        thisData = theseScores*allCoeff{c}';
        thisData = mean(thisData,1);
        [~,thisMin] = min(thisData);
        plot(thisData,'Color',plotColours(w,:),'LineWidth',plotLineWidth);
        xline(thisMin,'Color',plotColours(w,:),'LineWidth',plotLineWidth);
    end
end

%% Manuscript Figure
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.10 0.08], [0.1 0.06], [0.08 0.02]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis.
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins
%

poolParticipants = 0;

makefigure(17.17,9);

plot_settings;
whichSlots = [1 2 3];

for c = 1:3
    subplot(2,3,whichSlots(c))
    
    if poolParticipants
        % Plot coeff from pooled data
        thisPlot = plot(trialTimes{c},allCoeff{c}(:,1:2),'LineWidth',plotLineWidth);
        for pcaI = 1:2
            thisPlot(pcaI).Color = plotColours(pcaI,:);
        end
        
    else
        % Plot mean coeff (unpooled)
        smoothPCACoeff{c} = movmean(allPCACoeff{c},mmeanwindow,2);
        meanPCACoeff = squeeze(mean(smoothPCACoeff{c},1));
        stdPCACoeff = squeeze(std(smoothPCACoeff{c},[],1));
        ciPCACoeff = tVal*stdPCACoeff./sqrt(nParticipants);
        b1 = boundedline(trialTimes{c},meanPCACoeff(:,1),ciPCACoeff(:,1),'cmap',plotColours(1,:), 'alpha');
        b1.LineWidth = plotLineWidth;
        hold on;
        b2 = boundedline(trialTimes{c},meanPCACoeff(:,2),ciPCACoeff(:,2),'cmap',plotColours(2,:), 'alpha');
        b2.LineWidth = plotLineWidth;
    end
    
    ax = gca;
    set_axis(ax, axisSettings);
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'Comp. amp. (a.u.)';
    ax.XLim = [0 2.5];
    if poolParticipants
        l1 = legend('PC1','PC2');
    else
        l1 = legend([b1,b2],'PC1','PC2');
    end
    l1.Box = 'off';
    l1.Location = 'best';
    
    if ismember(c,[1 2 3])
        ax.XLabel.String = '';
    end
    
    if ismember(c,[2 3])
        ax.YLabel.String = '';
    end
end

meanScores = squeeze(mean(allPCAScores(:,:,whichChannel,:),2));

% Effect of PCA2
whichPlots = [4 5 6];
plot_settings;
pc2ScoresToPlot = [];
allMinMax = [];
for c = 1:3
    if poolParticipants
        thisMin = -25;
        thisMax = 25;
    else
        thisMin = round(mean(min(allPCAScores2{c}(:,:,2),[],2)));
        thisMax = round(mean(max(allPCAScores2{c}(:,:,2),[],2)));
        allMinMax(c,1) = thisMin;
        allMinMax(c,2) = thisMax;
    end
    pc2ScoresToPlot = [pc2ScoresToPlot; linspace(thisMin,thisMax,7)];
end
disp('Min/Max PC2 scores for the figure');
disp(allMinMax);

for c = 1:3
    subplot(2,3,whichPlots(c));
    plotColours = cbrewer('div','RdYlBu',size(pc2ScoresToPlot,2));
    plotColours = cbrewer('div','RdBu',size(pc2ScoresToPlot,2));
    plotColours = flip(plotColours);
    
    plotHandles = [];
    meanPC1Shape = mean(theseScores(:,1)*allCoeff{c}(:,1)',1);
    
    % Indiviual PCAs
    meanPCACoeffInd = squeeze(mean(smoothPCACoeff{c},1));
    stdPCACoeffInd = squeeze(std(smoothPCACoeff{c},[],1));
    ciPCACoeffInd = tVal*stdPCACoeffInd./sqrt(nParticipants);
    
    meanPC1ShapeInd = [];
    meanPC2CoeffInd = [];
    for iParticipant = 1:nParticipants
        thisPC1Coeff = smoothPCACoeff{c}(iParticipant,:,1);
        thisPC1Score = squeeze(allPCAScores2{c}(iParticipant,:,1));
        thisPC1Rec = thisPC1Score' * thisPC1Coeff;
        thisPC1Rec = mean(thisPC1Rec,1);
        meanPC1ShapeInd = [meanPC1ShapeInd; thisPC1Rec];
        
        thisPC2Coeff = smoothPCACoeff{c}(iParticipant,:,2);
        meanPC2CoeffInd = [meanPC2CoeffInd; thisPC2Coeff];
    end
    grandMeanPC1ShapeInd = mean(meanPC1ShapeInd,1);
    grandMeanPC2CoeffInd = mean(meanPC2CoeffInd,1);
    
    for w = 1:size(pc2ScoresToPlot,2)
        if poolParticipants
            thisPC2Shape = mean(pc2ScoresToPlot(c,w)*allCoeff{c}(:,2)',1);
            thisData = meanPC1Shape + thisPC2Shape;
        else
            % Individual
            thisPC2Shape = pc2ScoresToPlot(c,w)*grandMeanPC2CoeffInd;
            thisData = grandMeanPC1ShapeInd + thisPC2Shape;
            
        end
        
        [~,thisMin] = min(thisData);
        plotHandles(w) = plot(trialTimes{c},thisData,'Color',plotColours(w,:),'LineWidth',plotLineWidth);
        hold on;
        xlabel('Time (s)');
        ylabel('Voltage (\muV)');
        ax = gca;
        set_axis(ax, axisSettings);
        if poolParticipants
            % ax.YLim(1) = ax.YLim(1) - 0.5;
            % ax.YLim = [-1 0.5];
        else
            ax.YLim  = [-15 2];
        end
        ax.XLim = [0 2.5];
    end
    if poolParticipants
        plotHandles(size(pc2ScoresToPlot,2)+1) = plot(trialTimes{c},meanPC1Shape,'k','LineWidth',plotLineWidth);
    else
        % Individual
        plotHandles(size(pc2ScoresToPlot,2)+1) = plot(trialTimes{c},grandMeanPC1ShapeInd,'k','LineWidth',plotLineWidth);
    end
    
    if ismember(c,[2 3])
        ax.YLabel.String = '';
    end
end

save(fullfile(resultsFolder,'pca_results_production_task.mat'));
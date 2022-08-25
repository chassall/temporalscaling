% Make rERPs and ERPs in the Prediction task
%
% Other m-files required: eeglab

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% 25-Aug-2022
% Requires pinv_reg, makeStretchBasis, find_artifacts, basicrap.m

eeglab; close all; clear all; rng(41); % for reproducibility

% Set data and results folders (change as needed)
dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task3_prediction/data';
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

ps = [1:12 14 16:21]; % 13, 15 missing on datadryad
participants = num2str(ps,'%0.2d,');
participants(end) = [];
participants = strsplit(participants,',');
nParticipants = length(participants);
srate = 200;
doLaplacian = 0;
task = 'entrainrhythm1b';

% Rhythmic and repeated
whichConditions = {{131,231},{142,242},{133,233},{144,244}}; % s cue, l cue, s target, l target

% Artifact settings
maxMin = 150;
level = 150;
step = 50;
lowest = 0.1;

warningInterval = [-0.2 0.8];
targetInterval = [-0.8 0.2];
warningTimes = warningInterval(1):1/srate:warningInterval(2)-1/srate;
targetTimes = targetInterval(1):1/srate:targetInterval(2)-1/srate;
interval = {warningInterval,warningInterval,targetInterval,targetInterval};
baseline = {[-200 0],[-200 0],[-20 20],[-20 20]};
numPoints = round(srate*diff(interval{1}));
erps = nan(length(ps),4,64,numPoints);

longCueInterval = [-0.2 1.3];
longCueTimes = longCueInterval(1):1/srate:longCueInterval(2)-1/srate;
numLongCuePoints = round(srate*diff(longCueInterval));
longCueERPs = nan(length(ps),64,numLongCuePoints);

longTargetInterval = [-1.3 0.2];
longTargetTimes = longTargetInterval(1):1/srate:longTargetInterval(2)-1/srate;
numLongTargetPoints = round(srate*diff(longTargetInterval));
longTargetERPs = nan(length(ps),64,numLongTargetPoints);

% GLM parameters
doReg = 1;
amplitudeThreshold = level; % Artifact threshold
winms = 2000; % Artifact window
stepms = 100; % Artifact detection step size
tsMethod = 1; %
scaleLength = 1*200; % 0.7*200, 1.3*200, Length of temporally-scaled (TS) component

testMargin = 0;
margin = testMargin;
warningMargin = testMargin;
targetMargin = testMargin;

stretchMethod = 'box';

allWarningBetas = nan(length(ps),64,numPoints);
allTargetBetas = nan(length(ps),64,numPoints);
allTSBetas = nan(length(ps),64,scaleLength);
allRSquared = nan(length(participants),64);
allRSquaredFixed = nan(length(participants),64);

allArtifacts = zeros(length(ps),3);
allArtifactsByChannel = zeros(length(ps),4,64);

% Mean proportion of artifacts of either type
allContinuousArtifacts = nan(length(ps),1);

% Tikhonov parameters
loadPrevReg = 1;
lambdas = [0.001 0.01 1 10 100 1000 10000 100000 1000000 10000000];
k = 10; % k-fold cross-validation
regChannels = {'FCz'}; % DO the cross-validation at these channels
tsBreakpoints = [numPoints (numPoints+numPoints)]; % Since the betas will be appended, we need to know where the boundaries are

notsBreakpoints = [200];
regtype = 'onediff';
if ~loadPrevReg
    allCVErrors = nan(length(ps),64,length(lambdas));
else
    load('predictionCVErrors.mat');
end

allX = {};

for iParticipant = 1:length(ps)

        preprocessedFolder = fullfile(dataFolder,'derivatives','eegprep',['sub-' participants{iParticipant}]);
        preprocessedFile = ['sub-' participants{iParticipant} '_task-' task '_eegprep'];
        load(fullfile(preprocessedFolder,preprocessedFile), 'EEG'); 
        eyeLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Eye'));
        [~,I] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
        whichOnes = find(I == eyeLabel); % Max out of all possible labels
        EEG = pop_subcomp(EEG,whichOnes,0);

    allEpochedEEG = {};
    for w = 1:length(whichConditions)
        allEpochedEEG{w} = pop_epoch(EEG,whichConditions{w},interval{w});
        allEpochedEEG{w} = pop_rmbase(allEpochedEEG{w},baseline{w});
        [theseArtifacts,isArtifactsCT] = find_artifacts(allEpochedEEG{w}, maxMin, level, step, lowest);
        channelPercent = squeeze(mean(isArtifactsCT,2));
        allArtifactsByChannel(iParticipant,w,:) = channelPercent;
        thisERP = mean(allEpochedEEG{w}.data(:,:,~theseArtifacts),3);
        erps(iParticipant,w,:,:) = thisERP;
    end
    
    wEEG = pop_epoch(EEG,[whichConditions{1:2}],interval{1});
    tEEG = pop_epoch(EEG,[whichConditions{3:4}],interval{3});
    sEEG = pop_epoch(EEG,whichConditions{1},interval{1});
    
    wArtifacts = find_artifacts(wEEG, maxMin, level, step, lowest);
    tArtifacts = find_artifacts(tEEG, maxMin, level, step, lowest);
    sArtifacts = find_artifacts(sEEG, maxMin, level, step, lowest);
    
    longCueEEG = pop_epoch(EEG,whichConditions{2},longCueInterval);
    longCueEEG = pop_rmbase(longCueEEG,baseline{2});
    longCueArtifacts = find_artifacts(longCueEEG, maxMin, level, step, lowest);
    longCueERPs(iParticipant,:,:) = mean(longCueEEG.data(:,:,~longCueArtifacts),3);
    longCueEEG = pop_epoch(EEG,whichConditions{2},longCueInterval);
    
    longTargetEEG = pop_epoch(EEG,whichConditions{4},longTargetInterval);
    longTargetEEG = pop_rmbase(longTargetEEG,baseline{4});
    longTargetArtifacts = find_artifacts(longTargetEEG, maxMin, level, step, lowest);
    longTargetERPs(iParticipant,:,:) = mean(longTargetEEG.data(:,:,~longTargetArtifacts),3);
    
    if doReg
        
        % Artifact detection
        [WinRej, chanrej] = basicrap(EEG, 1:EEG.nbchan, amplitudeThreshold, winms, stepms, 1,[],[],'peak-to-peak',1);
        
        %isOK = ~any(thisLatency >= WinRej(:,1) & thisLatency <= WinRej(:,2));
        eegTypes = [EEG.event.type];
        
        % Indices
        short_warning_i = find(ismember(eegTypes,[whichConditions{1}{:}]));
        short_warning_times = round([EEG.event(short_warning_i).latency]);
        
        long_warning_i = find(ismember(eegTypes,[whichConditions{2}{:}]));
        long_warning_times = round([EEG.event(long_warning_i).latency]);
        
        warning_i = find(ismember([EEG.event.type],[whichConditions{1}{:} whichConditions{2}{:}]));
        target_i = find(ismember([EEG.event.type],[whichConditions{3}{:} whichConditions{4}{:}]));
        
        warning_times = round([EEG.event(warning_i).latency]);
        target_times = round([EEG.event(target_i).latency]);
        
        % Regressors
        r1 = sparse(EEG.pnts,numPoints); % Warning
        r2 = sparse(EEG.pnts,numPoints); % Target
        r3 = sparse(EEG.pnts,scaleLength); % TS
        
        for w = 1:length(warning_times)
            r1(warning_times(w)+srate*interval{1}(1):warning_times(w)+srate*interval{1}(2)-1,:) = makeStretchedBasis(numPoints,numPoints,1,warningMargin,'box');
        end
        
        for t = 1:length(target_times)
            r2(target_times(t)+srate*interval{3}(1):target_times(t)+srate*interval{3}(2)-1,:) = makeStretchedBasis(numPoints,numPoints,1,targetMargin,'box');
        end
        
        %% regressor r3 = time-scaled

        for t = 1:length(short_warning_times)
            thisBeepTime = short_warning_times(t);
            thisResponseTime = thisBeepTime + 0.7*EEG.srate;
            r3(thisBeepTime:thisResponseTime,:) = makeStretchedBasis(thisResponseTime-thisBeepTime+1,scaleLength,tsMethod,margin,stretchMethod);
        end
        
        for t = 1:length(long_warning_times)
            thisBeepTime = long_warning_times(t);
            thisResponseTime = thisBeepTime + 1.3*EEG.srate;
            r3(thisBeepTime:thisResponseTime,:) = makeStretchedBasis(thisResponseTime-thisBeepTime+1,scaleLength,tsMethod,margin,stretchMethod);
        end
        
        %% Make design matrix
        X = [r1 r2 r3];
        XFixed = [r1 r2];
        hasFixed = any(XFixed,2);

        XNoTS = [r1 r2]; % Fixed-time components only (no time-scaling)
        XNoFixed = [r3];
        if doLaplacian
            [data,G,H] = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
        else
            data = EEG.data;
        end
        
        % Store an copy for the purpose of constucting the residuals
        oX = X; % Original X
        oData = data; % Original EEG
        
        allX{iParticipant} = X;
        
        % Get Artifact indices
        iArtifact = [];
        for a = 1:size(WinRej,1)
            iArtifact = [iArtifact WinRej(a,1):WinRej(a,2)];
        end
        toExclude = iArtifact;
        
        % Determine empty rows
        isEmptyRow = sum(abs(X),2) == 0;
        
        % Determine artifact row
        isArtifactRow = false(size(isEmptyRow));
        isArtifactRow(iArtifact) = true;
        
        % Proportion of included rows that are "bad"
        isIncludedArtifact = ~isEmptyRow & isArtifactRow;
        allContinuousArtifacts(iParticipant) = mean(full(isIncludedArtifact));
        
        % Now actually remove empty/artifact rows
        toRemove = isEmptyRow | isArtifactRow;
        X(toRemove,:) = [];
        data(:,toRemove) = [];
        XFixed(toRemove,:) = [];
        hasFixed(toRemove) = [];

        iChannel = eeg_chaninds(EEG,regChannels);
        % Try Tikhonov regression
        if ~loadPrevReg
            allErrors = nan(EEG.nbchan,length(lambdas));
            tempBetas = [];
            for li = 1:length(lambdas)
                lambda = lambdas(li);
                cv = cvpartition(size(data,2),'KFold',k);
                theseErrors = nan(1,k);
                
                for c = iChannel
                    for ki = 1:k
                        thisTrainingI = training(cv,ki);
                        thisTestingI = test(cv,ki);
                        thisPinv = pinv_reg(X(thisTrainingI,:),lambda,regtype, tsBreakpoints);
                        thisBeta = thisPinv * data(c,thisTrainingI)';
                        
                        % For testing
                        thisResidual = data(c,thisTestingI)' - full(X(thisTestingI,:)) * thisBeta;
                        thisError = sum(thisResidual.*thisResidual);
                        theseErrors(ki) = thisError;
                    end
                    
                    allErrors(c,li) = mean(theseErrors);
                end
            end
            allCVErrors(iParticipant,:,:) = allErrors;
            
            plot(allErrors' - mean(allErrors'));
            title(num2str(iParticipant));
            drawnow();
        else
            allErrors = squeeze(allCVErrors(iParticipant,:,:));
        end
        
        % Examine errors to pick the best lambda
        [M,I] = min(allErrors');
        bestLambdaI = mode(I);
        bestLambda = lambdas(I(iChannel));
        
        % Do regression
        pDM = pinv_reg(X,bestLambda,regtype, tsBreakpoints);
        pDMFixed = pinv_reg(XFixed,bestLambda,regtype,tsBreakpoints(1));
        for c = 1:EEG.nbchan
            thisChannelData =  double(data(c,:));
            thisBeta = pDM * thisChannelData';
            beta(:,c) = thisBeta;
            allBetasStim(iParticipant,c,:) = thisBeta(1:numPoints);
            allBetasResp(iParticipant,c,:) = thisBeta(numPoints+1:2*numPoints);
            allBetasTS(iParticipant,c,:) = thisBeta(2*numPoints+1:2*numPoints+scaleLength);
            betaFixed(:,c) = pDMFixed * thisChannelData';
        end

        % Compute error
        residual = data - (X * beta)';
        residual(:,~hasFixed) = []; % Only count samples that are represented in both analyses (to make it fair)
        rsquared = sum(residual.* residual,2);
        residualFixed = data - (XFixed*betaFixed)';
        residualFixed(:,~hasFixed) = [];
        rquaredFixed = sum(residualFixed.* residualFixed,2);
        allRSquared(iParticipant,:) = rsquared;
        allRSquaredFixed(iParticipant,:) = rquaredFixed;

        % Compute residual datasets
        residualsWithTS = oData - (XNoTS * beta(1:size(XNoTS,2),:))';
        rEEGWithTS.data = residualsWithTS; % TS only
        residualsWithFixed = oData - (XNoFixed * beta(size(XNoTS,2)+1:end,:))';
        rEEGWithFixed.data = residualsWithFixed; % Fixed only
        
    end
    
    % Save the residuals for this participant
    if doReg || loadPrevReg
        if doLaplacian
            residualsFolder = fullfile(dataFolder,'derivatives','eegresiduallap',['sub-' participants{iParticipant}]);
        else
            residualsFolder = fullfile(dataFolder,'derivatives','eegresidual',['sub-' participants{iParticipant}]);
        end
        if ~exist(residualsFolder,'dir'), mkdir(residualsFolder); end
        residualsFile = ['sub-' participants{iParticipant} '_task-prediction_eegresidual' stretchMethod];
        save(fullfile(residualsFolder,residualsFile),'rEEGWithTS','rEEGWithFixed','beta','toExclude');
    end
    
    if ~doReg
        if doLaplacian
            residualsFolder = fullfile(dataFolder,'derivatives','eegresiduallap',['sub-' participants{iParticipant}]);
        else
            residualsFolder = fullfile(dataFolder,'derivatives','eegresidual',['sub-' participants{iParticipant}]);
        end
        if ~exist(residualsFolder,'dir'), mkdir(residualsFolder); end
        residualsFile = ['sub-' participants{iParticipant} '_task-prediction_eegresidual_noreg'];
        save(fullfile(residualsFolder,residualsFile),'rEEGWithTS','rEEGWithFixed','beta','toExclude');
    end
end

close all;
save(fullfile(resultsFolder,['rerp_and_erp_results_prediction_task_' num2str(doLaplacian)  '.mat']));

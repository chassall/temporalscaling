% Plots the rERPs
%
% Other m-files required: subtightplot, cbrewer, eeglab, boundedline

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% June 2020; Last revision: 18-Jun-2020

close all; clear all; clc; rng(41); % for reproducibility

% resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';
resultsFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\results';
tiffRes = '-r600';

whichTask = 'perception';
whichInterpMethod = 'box';
tiffFolder = './tiffs/';
doLaplacian = 0;
switch whichTask
    case 'production'
        numParticipants = 20;
        load(fullfile(resultsFolder,['rerp_results_production_task_' whichInterpMethod '_' num2str(doLaplacian) '.mat']),'respBL','stimBL','respScalePnts','stimScalePnts','scaleLength','lambdas','allCVErrors','participants','allContinuousArtifacts','allBetasStim','allBetasResp','allBetasTS','EEG','allRSquared','allRSquaredFixed');
    case 'perception'
        numParticipants = 20;
        load(fullfile(resultsFolder,['rerp_results_perception_task_'  whichInterpMethod '_' num2str(doLaplacian)  '.mat']),'respBL','stimBL','respScalePnts','stimScalePnts','scaleLength','lambdas','allCVErrors','participants','allContinuousArtifacts','allBetasStim','allBetasResp','allBetasTS','EEG','allRSquared','allRSquaredFixed');
    case 'prediction'
        numParticipants = 19;
        load(fullfile(resultsFolder,['rerp_and_erp_results_prediction_task_' num2str(doLaplacian) '.mat']),'targetInterval','srate','warningInterval','ps','respBL','stimBL','respScalePnts','stimScalePnts','scaleLength','lambdas','allCVErrors','participants','allContinuousArtifacts','allBetasStim','allBetasResp','allBetasTS','EEG','allRSquared','allRSquaredFixed');
    case 'dm'
        numParticipants = 18;
        load(fullfile(resultsFolder,['rerp_results_dm_task_reref_' num2str(doLaplacian) '.mat']));
        % load('../results/rerp_results_dm_task_0.mat');
end
plot_settings;

%% Artifacts
disp('Artifacts');
meanA = mean(allContinuousArtifacts);
stdA = std(allContinuousArtifacts);
tVal = abs(tinv(0.025,length(participants)-1));
ciA = tVal * stdA / length(participants);
disp([meanA meanA-ciA meanA+ciA]);

%% Plot best lambdas across all channels

subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.06], [0.14 0.1], [0.04 0.01]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

makefigure(48,20);
for iParticipant = 1:numParticipants
   subplot(4,5,iParticipant);
   
   allErrors = squeeze(allCVErrors(iParticipant,:,:));
   plot(allErrors' - mean(allErrors'),'LineWidth',plotLineWidth);
   xlabel('\lambda');
   ylabel('SSR');
   ax = gca;
   ax.Box = 'off';
   ax.XTickLabel = lambdas;
   % ax.XLim = [1 6];
end
switch whichTask
    case 'production'
        theseErrors = squeeze(allCVErrors(:,10,:));
        [M,I] = min(theseErrors,[],2);
    case 'perception'
        theseErrors = squeeze(allCVErrors(:,10,:));
        [M,I] = min(theseErrors,[],2);
    case 'prediction'
        theseErrors = squeeze(allCVErrors(:,47,:));
        [M,I] = min(theseErrors,[],2);
    case 'dm'
        theseErrors = squeeze(allCVErrors(:,19,:));
        [M,I] = min(theseErrors,[],2);
end

disp('Best Lambdas');
disp(lambdas(I));

%% R-Squared (Fixed-only verus Scaled + Fixed)
% Average across all channels
meanError = mean(allRSquared(:,:),2);
meanErrorFixed = mean(allRSquaredFixed(:,:),2);

disp('MSE, Full Model');
mean(meanError)
[H,P,CI,STATS] = ttest(meanError);
CI

disp('MSE, Fixed-Onlu Model');
mean(meanErrorFixed)
[H,P,CI,STATS] = ttest(meanErrorFixed);
CI

disp('R-Squared, Fixed+Scaled Versus Fixed-Only')
meanScores = meanError - meanErrorFixed;
sdScores = std(meanScores);
[H,P,CI,STATS] = ttest(meanScores);
disp(['Mean: ' num2str(mean(meanScores))]);
disp(['CI: ' num2str(CI')]);
cohensD = mean(meanScores)/STATS.sd;
disp('R-Squared Component Stats (t,p,Cohen''s d)');
disp([STATS.tstat,P,cohensD]);
disp('Shapiro-Wilk test');
[H, pValue, SWstatistic] = swtest(meanScores)

%% Manuscript Figure 2

switch whichTask
    case 'production'
        stimChannelString = 'FCz';
        respChannelString = 'FCz';
        tsChannelString =  'P3';
        stimWindow = [0.084 0.124];
        respWindow = [0 0.040];
        alphaValue = 0.001;
    case 'perception'
        stimChannelString = 'F4';
        respChannelString = 'FCz';
        tsChannelString =  'Cz';
        tsChannelString =  'FC1';
        stimWindow = [0.084 0.124]; % N1
        respWindow = [0.084 0.124];
         alphaValue = 0.01;
    case 'prediction'
        stimChannelString = 'P2';
        respChannelString = 'CP1';
        tsChannelString =  'Cz';
        stimWindow = [0.350 0.550];
        respWindow = [-0.06 0.135];
         alphaValue = 0.05;
        
        % "Predict" script output is a bit different, so...
        participants = ps;
        
        % Stim-locked component
        stimTimelimits = warningInterval; % Time window
        stimScalePntLimits = srate * stimTimelimits;
        stimScalePnts = srate * stimTimelimits(1):srate *(stimTimelimits(2)-1/srate); % TODO: make consistent
        numStimPnts = length(stimScalePnts);
        stimBL = [-0.2 0];
        
        % Response-locked component
        respTimelimits = targetInterval; % Time window
        respScalePntLimits = srate * respTimelimits;
        respScalePnts = srate * respTimelimits(1): srate *(respTimelimits(2)-1/srate); % TODO: make consistent
        numRespPnts = length(respScalePnts);
        respBL = [-0.8, -0.6];
        
        % Time-scaled component
        scaleLength = 200; % In data points
        margin = 0; % Num points to exclude from endpoints
    case 'dm'
        stimChannelString = 'Pz';
        respChannelString = 'FC1';
        cppChannelString = 'Pz';
        tsChannelString =  'Cz';
        stimWindow = [0.455 0.610]; % N1
        respWindow = [-0.025 0.065];
         alphaValue = 0.05;
end


tVal = abs(tinv(0.025,length(participants)-1));
stimTimes = stimScalePnts/EEG.srate;
respTimes = respScalePnts/EEG.srate;
tsBL = scaleLength*.05; % Size of baseline, in samples (baseline at end of interval for now)

switch whichTask
    case {'production'}
        yLimits = [-12 5];
    case {'perception'}
        yLimits = [-12 5];
    case 'prediction'
        yLimits = [-16 25];
    case 'dm'
        yLimits = [-3 6.5];
end

stimBLI = dsearchn(stimTimes',stimBL');
respBLI = dsearchn(respTimes',respBL');

allBetasStim(:,:,:) = allBetasStim(:,:,:) - nanmean(allBetasStim(:,:,stimBLI(1):stimBLI(2)),3);
allBetasResp(:,:,:) = allBetasResp(:,:,:) - nanmean(allBetasResp(:,:,respBLI(1):respBLI(2)),3);
allBetasTS(:,:,:) = allBetasTS(:,:,:) - nanmean(allBetasTS(:,:,1:tsBL),3);

stimChannel = eeg_chaninds(EEG,{stimChannelString});
stim = squeeze(nanmean(allBetasStim(:,stimChannel,:),1));
stimCI = squeeze(nanstd(allBetasStim(:,stimChannel,:),[],1));
stimCI = tVal * stimCI ./ sqrt(length(participants));

respChannel = eeg_chaninds(EEG,{respChannelString});
resp = squeeze(nanmean(allBetasResp(:,respChannel,:),1));
respCI = squeeze(nanstd(allBetasResp(:,respChannel,:),[],1));
respCI = tVal * respCI ./ sqrt(length(participants));

if strcmp(whichTask,'dm')
cppChannel = eeg_chaninds(EEG,{cppChannelString});
cpp = squeeze(nanmean(allBetasResp(:,cppChannel,:),1));
cppCI = squeeze(nanstd(allBetasResp(:,cppChannel,:),[],1));
cppCI = tVal * cppCI ./ sqrt(length(participants));
end

tsChannel = eeg_chaninds(EEG,{tsChannelString});
ts = squeeze(nanmean(allBetasTS(:,tsChannel,:),1));
tsCI = squeeze(nanstd(allBetasTS(:,tsChannel,:),[],1));
tsCI = tVal * tsCI ./ sqrt(length(participants));

%% Permutation tests (time, channel) for revision 2

% Load fieldtrip template for defining neighbours and delete the ones that
% we don't need for this task
load('elec1010_neighb.mat');
ourChannels = {EEG.chanlocs.labels};
toDelete = [];
for i = 1:length(neighbours)
    invNeighb = ~ismember(neighbours(i).neighblabel,ourChannels);
    neighbours(i).neighblabel(invNeighb) = [];
    thisLabel = neighbours(i).label;
    if ~ismember(thisLabel,ourChannels)
        toDelete = [toDelete i];
    end
end
neighbours(toDelete) = []; % Remove invalid channels

% Compute actual t values for each channel and time
[actualhs, ps, cis, stats] = ttest(allBetasTS,zeros(size(allBetasTS)),'Alpha',alphaValue);
actualhs = squeeze(actualhs);
actualts = squeeze(stats.tstat);

% Now loop through each channel and compute a cluster statistic based on
% the temporal and spatial adjacency
clusterInfo = {};
for i = 1:length(neighbours)

    % This channel
    thisLabel = neighbours(i).label;
    clusterInfo{i,1} = thisLabel;
    thisI = eeg_chaninds(EEG,thisLabel);
    thisH = actualhs(thisI,:);
    thisT = actualts(thisI,:);

    % Neighbouring channels
    theseNeighbours = neighbours(i).neighblabel;
    neighbI = eeg_chaninds(EEG,theseNeighbours);
    neighbH = actualhs(neighbI,:);
    neighbT = actualts(neighbI,:);

    allH = [thisH; neighbH];
    commonH = all(allH,1);

    allT = [thisT; neighbT];

    % Find temporal clusters
    isDiff = [true; diff(commonH(:)) ~= 0];
    commonH = [commonH 0];
    changeInd = find([isDiff' true]);
    rejectNull = commonH(changeInd);
    actualClusterExtents = diff(changeInd);
    actualClusterExtents = actualClusterExtents(rejectNull(1:end-1) == 1); % Only count clusters where h = 1
    
    actualclusterStartI = changeInd(rejectNull == 1);
    actualclusterEndI = actualclusterStartI + actualClusterExtents - 1;
    actualclusterMasses = [];
    for j = 1:length(actualClusterExtents)
        actualclusterMasses(j) = mean(sum(allT(:,actualclusterStartI(j):actualclusterEndI(j)),2),1);
    end
    [actualMaxCluster,maxI] = max(abs(actualclusterMasses));
    scaledWindow = [actualclusterStartI(maxI) actualclusterEndI(maxI)];

    clusterInfo{i,2} = actualMaxCluster;
    clusterInfo{i,3} = scaledWindow;
    clusterInfo{i,4} = scaledWindow./size(allBetasTS,3);
end

%% Permutation tests (new)

numPerms = 1000;
allPerms = -1 + 2*(randi(2,size(allBetasTS,1),numPerms)-1);
clusterMasses = [];
for i = 1:size(allPerms,2)
    thisPerm = allPerms(:,i);
    flippedBetas = allBetasTS .* thisPerm;
    [hs, ps, cis, stats] = ttest(flippedBetas,zeros(size(flippedBetas)),'Alpha',alphaValue);
    
    hs = squeeze(hs);
    tstats = squeeze(stats.tstat);

    for i = 1:length(neighbours)

        % This channel
        thisLabel = neighbours(i).label;
        thisI = eeg_chaninds(EEG,thisLabel);
        thisH = hs(thisI,:);
        thisT = actualts(thisI,:);

        % Neighbouring channels
        theseNeighbours = neighbours(i).neighblabel;
        neighbI = eeg_chaninds(EEG,theseNeighbours);
        neighbH = hs(neighbI,:);

        allH = [thisH; neighbH];
        commonH = all(allH,1);

        allT = [thisT; neighbT];

        % Find temporal clusters
        isDiff = [true; diff(commonH(:)) ~= 0];
        commonH = [commonH 0];
        changeInd = find([isDiff' true]);
        rejectNull = commonH(changeInd);
        actualClusterExtents = diff(changeInd);
        actualClusterExtents = actualClusterExtents(rejectNull(1:end-1) == 1); % Only count clusters where h = 1

        actualclusterStartI = changeInd(rejectNull == 1);
        actualclusterEndI = actualclusterStartI + actualClusterExtents - 1;
        theseClusterMasses = [];
        for j = 1:length(actualClusterExtents)
            theseClusterMasses(j) = mean(sum(allT(:,actualclusterStartI(j):actualclusterEndI(j)),2),1);
        end
        [actualMaxCluster,maxI] = max(abs(theseClusterMasses));

        if isempty(actualMaxCluster)
            clusterMasses = [clusterMasses 0];
        else
            clusterMasses = [clusterMasses actualMaxCluster];
        end
    end

end

% Look at the distribution of clustermasses (unused)
figure();
histogram(clusterMasses);
for i = 1:size(clusterInfo,1)
   if ~isempty(clusterInfo{i,2})
   xline( clusterInfo{i,2}); hold on;
   end
end

% Check our clusters against the permutation
for i = 1:size(clusterInfo,1)
    if ~isempty(clusterInfo{i,2})
        clusterInfo{i,5} = sum(clusterMasses < clusterInfo{i,2}) / length(clusterMasses);
        clusterInfo{i,6} = 1 - clusterInfo{i,5};
    end 
end

% Display p value of the target clusgter.
% Can't use the index from the EEG file because it might be
% different from what's in clusterInfo.
thisTSIndex = find(strcmp(tsChannelString,clusterInfo(:,1))); 
disp('Cluster Info');
disp(clusterInfo(thisTSIndex,:));

disp('Neighbours');
disp(neighbours(thisTSIndex).label);
disp(neighbours(thisTSIndex).neighblabel);
clusterIs = [eeg_chaninds(EEG,neighbours(thisTSIndex).label) eeg_chaninds(EEG,neighbours(thisTSIndex).neighblabel)];

% %% Permutation tests (old, unused)
% 
% % Actual t-stats
% %tsBetaPByTime = squeeze(allBetasStim(:,tsChannel,:));
% %tsBetaPByTime = squeeze(allBetasResp(:,tsChannel,:));
% tsBetaPByTime = squeeze(allBetasTS(:,tsChannel,:));
% alphaValue = 0.05;
% % tCritical = tinv(alphaValue/2,size(allBetasTS,1)-1);
% [actualhs, ps, cis, stats] = ttest(tsBetaPByTime,zeros(size(tsBetaPByTime)),'Alpha',alphaValue);
% tstats = stats.tstat;
% isSig = find(actualhs);
% % isSig = abs(tstats) > abs(tCritical);
% isDiff = [true; diff(actualhs(:)) ~= 0];
% actualhs = [actualhs 0];
% changeInd = find([isDiff' true]);
% rejectNull = actualhs(changeInd);
% actualClusterExtents = diff(changeInd);
% actualClusterExtents = actualClusterExtents(rejectNull(1:end-1) == 1); % Only count clusters where h = 1
% 
% actualclusterStartI = changeInd(rejectNull == 1);
% actualclusterEndI = actualclusterStartI + actualClusterExtents - 1;
% actualclusterMasses = [];
% actualClusterHeights = [];
% for i = 1:length(actualClusterExtents)
%     actualclusterMasses(i) = sum(tstats(actualclusterStartI(i):actualclusterEndI(i)));
%     actualClusterHeights(i) = max(abs(tstats(actualclusterStartI(i):actualclusterEndI(i))));
% end
% [actualMaxCluser,maxI] = max(abs(actualclusterMasses));
% scaledWindow = [actualclusterStartI(maxI) actualclusterEndI(maxI)];
% 
% [actualMaxClusterHeight,maxI] = max(actualClusterHeights);
% 
% subplot(1,2,1); plot(ts); hold on; plot(isSig,zeros(1,length(isSig)),'.');
% 
% numParticipants = size(tsBetaPByTime,1);
% 
% % Permutations
% % switch whichTask
% %     case {'production','perception'}
% %        allPerms = dec2bin(0:2^numParticipants-1)-'0'; % Makes all binary combinations
% %        allPerms(allPerms == 0) = -1;
% %        allPerms = allPerms';
% %     otherwise
%         numPerms = 1000;
%         allPerms = -1 + 2*(randi(2,size(tsBetaPByTime,1),numPerms)-1);
% % end
% 
% permClusterExtents = [];
% permClusterMasses = [];
% permClusterHeights = [];
% for i = 1:size(allPerms,2)
%     thisPerm = allPerms(:,i);
%     flippedBetas = tsBetaPByTime .* thisPerm;
%     [hs, ps, cis, stats] = ttest(flippedBetas,zeros(size(flippedBetas)),'Alpha',alphaValue);
%     tstats = stats.tstat;
%     
%     isDiff = [true; diff(hs(:)) ~= 0];
%     hs = [hs 0];
%     changeInd = find([isDiff' true]);
%     rejectNull = hs(changeInd);
%     runLengths = diff(changeInd);
%     runLengths = runLengths(rejectNull(1:end-1) == 1); % Only count clusters where h = 1
%     permClusterExtents = [permClusterExtents runLengths];
%     
%     clusterStartI = changeInd(rejectNull == 1);
%     clusterEndI = clusterStartI + runLengths - 1;
%     clusterMasses = [];
%     clusterHeights = [];
%     for i = 1:length(runLengths)
%         clusterMasses(i) = sum(tstats(clusterStartI(i):clusterEndI(i)));
%         clusterHeights(i) = max(abs(tstats(clusterStartI(i):clusterEndI(i))));
%     end
%     
%     [~,maxI] = max(abs(clusterMasses));
%     if isempty(clusterMasses)
%         permClusterMasses = [permClusterMasses 0];
%     else
%         permClusterMasses = [permClusterMasses clusterMasses(maxI)];
%     end
%     
%     [~,maxI] = max(clusterHeights);
%     permClusterHeights = [permClusterHeights clusterHeights(maxI)];
% end
% 
% subplot(1,2,2);
% histogram(abs(permClusterMasses));
% for i = 1:length(actualclusterMasses)
%    xline( abs(actualclusterMasses(i))); hold on;
% end
% 
% % Check our clusters against the permutation
% p = [];
% for i = 1:length(actualclusterMasses)
%     p(i) = sum(abs(permClusterMasses) < abs(actualclusterMasses(i))) / length(permClusterMasses);
% end
% disp(['Cluster ts: ' num2str(abs(actualclusterMasses))]);
% disp(['Cluster ps: ' num2str(1-p)]);
% 
% % p = [];
% % for i = 1:length(actualClusterExtents)
% %     p(i) = sum(permClusterExtents < actualClusterExtents(i)) / length(permClusterExtents);
% % end
% % disp(['Cluster lengths: ' num2str(actualClusterExtents)]);
% % disp(['Cluster ps: ' num2str(1-p)]);
% 
% % propLess = [];
% % for i = 1:length(actualClusterHeights)
% %     propLess(i) = sum(permClusterHeights < actualClusterHeights(i)) / length(permClusterExtents);
% % end
% % disp(['Cluster lengths: ' num2str(actualClusterHeights)]);
% % disp(['Cluster ps: ' num2str(1-propLess)]);
% % 
% % Greatest t cluster?
% disp(['(Points) From ' num2str(scaledWindow(1)) ' to ' num2str(scaledWindow(2))]);
% disp(['(Proportion) From ' num2str(scaledWindow(1)/size(allBetasTS,3)) ' to ' num2str(scaledWindow(2)/size(allBetasTS,3))]);


%% rERP Plot

% Set the analysis window for the scaled signal (from permutation testing
% above).
scaledWindow = clusterInfo{thisTSIndex,3};
scaledWindowProp = scaledWindow./length(ts);

disp('cluster extent');
disp(scaledWindowProp);

stimWinPnts = dsearchn(stimTimes',stimWindow');
respWinPnts = dsearchn(respTimes',respWindow');
scaledWinPnts = scaledWindow(1):scaledWindow(2);
areaAlpha = 0.5;

fig1 = makefigure(11,4);

subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.08], [0.2 0.15], [0.09 0.02]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

% Cue-Locked Plot
subplot(1,3,1);
yline(0); hold on;
bl = boundedline(stimTimes,stim,stimCI,'cmap',plotColours(1,:));
bl.LineWidth = plotLineWidth;
t1 = title('Cue'); t1.FontWeight = 'normal';
ax = gca;
ax.FontSize = plotFontSize;
ax.Box = 'off';
ax.YLabel.String = 'Voltage (\muV)';
ax.XLabel.String = 'Time (s)';
if ~doLaplacian
ax.YLim = yLimits;
end
ax.XLim = [stimTimes(1) stimTimes(end)];
text(ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))*0.05,ax.YLim(1)+ (ax.YLim(2)-ax.YLim(1))*0.1,stimChannelString,'FontSize',plotFontSize);
hold on;
%area(stimWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',areaAlpha,'LineStyle','none'); hold on;
%area(stimWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',areaAlpha,'LineStyle','none');
plot(stimWindow,[(ax.YLim(2) - ax.YLim(1))/10 (ax.YLim(2) - ax.YLim(1))/10],'Color',[0,0,0,areaAlpha],'LineWidth',3);

% Response/Probe/Target-Locked Plot
subplot(1,3,2);
yline(0);
bl = boundedline(respTimes,resp,respCI,'cmap',plotColours(2,:));
bl.LineWidth = plotLineWidth;
ax1 = gca;
switch whichTask
    case 'production'
        t2 = title('Response');
    case 'perception'
        t2 = title('Probe');
    case 'prediction'
        t2 = title('Target');
    case 'dm'
        t2 = title('Decision');
end
t2.FontWeight = 'normal';
ax = gca;
ax.FontSize = plotFontSize;
ax.Box = 'off';
ax.XLabel.String = 'Time (s)';
ax.XLim = [respTimes(1) respTimes(end)];
if ~doLaplacian
ax.YLim = yLimits;
end
hold on;
if strcmp(whichTask,'dm')
    b2 = boundedline(respTimes,cpp,cppCI,'cmap',plotColours(4,:));
    b2.LineWidth = plotLineWidth;
    ax2 = gca;
else
    text(ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))*0.05,ax.YLim(1)+ (ax.YLim(2)-ax.YLim(1))*0.1,respChannelString,'FontSize',plotFontSize);
end

%area(respWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',areaAlpha,'LineStyle','none'); hold on;
%area(respWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',areaAlpha,'LineStyle','none');
plot(respWindow,[(ax.YLim(2) - ax.YLim(1))/10 (ax.YLim(2) - ax.YLim(1))/10],'Color',[0,0,0,areaAlpha],'LineWidth',3);

if strcmp(whichTask,'dm')
    l = legend([b2,bl],cppChannelString,respChannelString,'box','off'); 
    l.Position = [0.4215    0.5228    0.1626    0.1803];
end
% Scaled Plot
subplot(1,3,3)
yline(0);
bl = boundedline(1:length(ts),ts,tsCI,'cmap',plotColours(3,:));
bl.LineWidth = plotLineWidth;
t3 = title('Scaled');
t3.FontWeight = 'normal';
ax = gca;
ax.FontSize = plotFontSize;
ax.Box = 'off';
ax.XLabel.String = 'Interval proportion (%)';
ax.XLim = [1 scaleLength];
ax.XTick = [1 scaleLength/2 scaleLength];
ax.XAxis.TickLabels = {'1','50','100'};
if ~doLaplacian
ax.YLim = yLimits;
end
text(ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))*0.05,ax.YLim(1)+ (ax.YLim(2)-ax.YLim(1))*0.1,tsChannelString,'FontSize',plotFontSize);
hold on;
%area(scaledWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',areaAlpha,'LineStyle','none'); hold on;
%area(scaledWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',areaAlpha,'LineStyle','none');
plot(scaledWindow,[(ax.YLim(2) - ax.YLim(1))/10 (ax.YLim(2) - ax.YLim(1))/10],'Color',[0,0,0,areaAlpha],'LineWidth',3);

%chi=get(gca, 'Children');
%set(gca, 'Children',[chi(end); chi(1:end-1)]);


% Print Figure
print([tiffFolder 'figure_3_' whichTask '_rerps_' whichInterpMethod '_' num2str(doLaplacian) '.tiff'],'-dtiff',tiffRes);

% Topos
stimTopo = mean(mean(allBetasStim(:,:,stimWinPnts(1):stimWinPnts(2)),3),1);
respTopo = mean(mean(allBetasResp(:,:,respWinPnts(1):respWinPnts(2)),3),1);
scaledTopo = mean(mean(allBetasTS(:,:,scaledWinPnts(1):scaledWinPnts(2)),3),1);
allTopo = [stimTopo; respTopo; scaledTopo];

% Stats
meanScores = mean(allBetasTS(:,tsChannel,scaledWinPnts(1):scaledWinPnts(2)),3);
sdScores = std(meanScores);
[H,P,CI,STATS] = ttest(meanScores);
disp(['Mean: ' num2str(mean(meanScores))]);
disp(['CI: ' num2str(CI')]);
cohensD = mean(meanScores)/STATS.sd;
disp('Scaled Component Stats (t,p,Cohen''s d)');
disp([STATS.tstat,P,cohensD]);
disp('Shapiro-Wilk test');
[H, pValue, SWstatistic] = swtest(meanScores)

%%

meanTS = squeeze(mean(allBetasTS,1));
minAcrossTime = min(meanTS,[],2);
[~,minChannelI] = min(minAcrossTime);
disp(['Where is the scaled signal minimum? ' EEG.chanlocs(minChannelI).labels]);

% Display max/min sensor locations
switch whichTask
    case 'production'
        [~,I] = min(stimTopo); EEG.chanlocs(I).labels
        [~,I] = min(respTopo); EEG.chanlocs(I).labels
        [~,I] = min(scaledTopo); EEG.chanlocs(I).labels
    case 'perception'
        [~,I] = min(stimTopo); EEG.chanlocs(I).labels
        [~,I] = min(respTopo); EEG.chanlocs(I).labels
        [~,I] = min(scaledTopo); EEG.chanlocs(I).labels
        
    case 'prediction'
        [~,I] = max(stimTopo); EEG.chanlocs(I).labels
        [~,I] = min(respTopo); EEG.chanlocs(I).labels
        [~,I] = min(scaledTopo); EEG.chanlocs(I).labels
    case 'dm'
        [~,I] = max(stimTopo); EEG.chanlocs(I).labels
        [~,I] = max(respTopo); EEG.chanlocs(I).labels
        [~,I] = min(scaledTopo); EEG.chanlocs(I).labels 
end


        
numContour = 5;
topoColour = cbrewer('div', 'RdBu', 256,'spline');
topoColour = flip(topoColour);
topoColour(topoColour < 0) = 0;
topoColour(topoColour > 1) = 1;
topoFigDim = [3,3];

% Set symmetrical topo limits
mapLimits = [floor(-max(max(abs(allTopo)))) ceil(max(max(abs(allTopo))))];

make_figure(topoFigDim(1),topoFigDim(2));
st = topoplot(stimTopo,EEG.chanlocs,'maplimits' ,mapLimits,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
%st = topoplot(stimTopo,EEG.chanlocs,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
st.Parent.XLim = [-0.6 0.6];
st.Parent.YLim = [-0.6 0.6];
set(gca,'color','none');
print([tiffFolder 'figure_3_' whichTask '_cue_topo_' whichInterpMethod '_' num2str(doLaplacian) '.tiff'],'-dtiff',tiffRes);

make_figure(topoFigDim(1),topoFigDim(2));
st = topoplot(respTopo,EEG.chanlocs,'maplimits' ,mapLimits,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
%st = topoplot(respTopo,EEG.chanlocs,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
st.Parent.XLim = [-0.6 0.6];
st.Parent.YLim = [-0.6 0.6];
print([tiffFolder 'figure_3_' whichTask '_response_topo_' whichInterpMethod '_' num2str(doLaplacian) '.tiff'],'-dtiff',tiffRes);

make_figure(topoFigDim(1),topoFigDim(2));
% st = topoplot([],EEG.chanlocs(clusterIs),'headcolor',[1 1 1],'style','blank'); hold on;
st = topoplot(scaledTopo,EEG.chanlocs,'maplimits' ,mapLimits,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
%st = topoplot(scaledTopo,EEG.chanlocs,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
st.Parent.XLim = [-0.6 0.6];
st.Parent.YLim = [-0.6 0.6];
print([tiffFolder 'figure_3_' whichTask '_scaled_topo_' whichInterpMethod '_' num2str(doLaplacian)  '.tiff'],'-dtiff',tiffRes);

make_figure(2.5,1);
c = colorbar();
colormap(topoColour);
c.TickLabels = [mapLimits(1) 0 mapLimits(2)];
c.Location = 'north';
c.Label.String = 'Voltage (\muV)';
c.Label.Visible = 'off';
c.FontSize = 6;
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
% ax.XAxisLocation = 'top';
print([tiffFolder 'figure_3_' whichTask '_colourbar_' whichInterpMethod '_' num2str(doLaplacian) '.tiff'],'-dtiff',tiffRes);

%% Manuscript Figures S4-6 (All Participant Waveforms)

participants = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};

axisSettings.fontSize = 8;
axisSettings.fontName = 'Arial';

subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.01 0.08], [0.04 0.02], [0.04 0.02]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins

make_figure(18.3,24.0);
nParticipants = numParticipants;
for iParticipant = 1:nParticipants
    subplot(nParticipants,3,(iParticipant-1)*3 + 1);
    plot(stimTimes, squeeze(allBetasStim(iParticipant,stimChannel,:)),'Color',plotColours(1,:),'LineWidth',plotLineWidth);
    ax = gca;
    set_axis(ax,axisSettings);
    if iParticipant == 1
        t1 = title(['Cue (' stimChannelString ')']);
    end
    ax.XLim = [stimTimes(1) stimTimes(end)];
    if iParticipant == nParticipants
        ax.YLabel.String = 'Voltage (\muV)';
        ax.XLabel.String = 'Time (s)';
    else
       ax.XAxis.Visible = 'off'; 
    end
    
    subplot(nParticipants,3,(iParticipant-1)*3 + 2);
    plot(respTimes,squeeze(allBetasResp(iParticipant,respChannel,:)),'Color',plotColours(2,:),'LineWidth',plotLineWidth);
    ax = gca;
    set_axis(ax,axisSettings);
    if iParticipant == 1
        switch whichTask
            case 'production'
                t2 = title(['Response (' respChannelString ')']);
            case 'perception'
                t2 = title(['Probe (' respChannelString ')']);
                
            case 'prediction'
                t2 = title(['Target (' respChannelString ')']);
            case 'dm'
                t2 = title(['Decision (' respChannelString ')']);
        end
    end
    ax.XLim = [respTimes(1) respTimes(end)];
    if iParticipant == nParticipants
        ax.YLabel.String = 'Voltage (\muV)';
        ax.XLabel.String = 'Time (s)';
    else
       ax.XAxis.Visible = 'off'; 
    end
    
    subplot(nParticipants,3,(iParticipant-1)*3 + 3);
    plot(squeeze(allBetasTS(iParticipant,tsChannel,:)),'Color',plotColours(3,:),'LineWidth',plotLineWidth);
    ax = gca;
    set_axis(ax,axisSettings);
    if iParticipant == 1
       title(['Scaled (' tsChannelString ')']); 
    end
    
    ax.XLim = [1 scaleLength];
    if iParticipant == nParticipants
        ax.YLabel.String = 'Voltage (\muV)';
        ax.XLabel.String = 'Interval proportion (%)';
        ax.XTick = [1 scaleLength/2 scaleLength];
        ax.XAxis.TickLabels = {'1','50','100'};
    else
        ax.XAxis.Visible = 'off';
    end
    
end

switch whichTask
    case 'production'
        print([tiffFolder 'sfig_5_' whichTask '_rerps_all_' whichInterpMethod '.tiff'],'-dtiff',tiffRes);
    case 'perception'
        print([tiffFolder 'sfig_6_' whichTask '_rerps_all.tiff'],'-dtiff',tiffRes);
    case 'prediction'
        print([tiffFolder 'sfig_7_' whichTask '_rerps_all.tiff'],'-dtiff',tiffRes);
    case 'dm'
        print([tiffFolder 'sfig_8_' whichTask '_rerps_all.tiff'],'-dtiff',tiffRes);
end

return;
%% Baseline correction/averaging

switch whichTask
    case 'production'
        stimChannelString = 'FCz';
        respChannelString = 'FCz';
        tsChannelString =  'Cz';
        stimWindow = [0.084 0.124]; % N1
        respWindow = [0 0.040];
        scaledWindow = round([2*scaleLength/5 4*scaleLength/5]);
        scaledWindow = round([6*scaleLength/10 8*scaleLength/10]);
    case 'perception'
        stimChannelString = 'FCz';
        respChannelString = 'FCz';
        tsChannelString =  'Cz';
        stimWindow = [0.084 0.124]; % N1
        respWindow = [0.084 0.124];
        scaledWindow = round([6*scaleLength/10 8*scaleLength/10]);
    case 'prediction'
        stimChannelString = 'POz';
        respChannelString = 'FCz';
        tsChannelString =  'CPz';
        stimWindow = [0.375 0.600];
        respWindow = [-0.360 -0.255];
        scaledWindow = round([6*scaleLength/10 8*scaleLength/10]);
        
        % "Predict" script output is a bit different, so...
        participants = ps;
        
        % Stim-locked component
        stimTimelimits = warningInterval; % Time window
        stimScalePntLimits = srate * stimTimelimits;
        stimScalePnts = srate * stimTimelimits(1):srate *(stimTimelimits(2)-1/srate); % TODO: make consistent
        numStimPnts = length(stimScalePnts);
        stimBL = [-0.2 0];
        
        % Response-locked component
        respTimelimits = targetInterval; % Time window
        respScalePntLimits = srate * respTimelimits;
        respScalePnts = srate * respTimelimits(1): srate *(respTimelimits(2)-1/srate); % TODO: make consistent
        numRespPnts = length(respScalePnts);
        respBL = [-0.02, 0.02];
        
        % Time-scaled component
        scaleLength = 200; % In data points
        margin = 0; % Num points to exclude from endpoints  
end

switch whichTask
    case {'production','perception'}
        yLimits = [-10 5];
    case 'prediction'
        yLimits = [-15 15];
end

switch whichTask
    
    case {'production','perception'}
        mapLimits = [-8 8];
    case 'prediction'
        mapLimits = [-12 12];
end

tVal = abs(tinv(0.025,length(participants)-1));

stimTimes = stimScalePnts/EEG.srate;
respTimes = respScalePnts/EEG.srate;

tsBL = 10; % Size of baseline, in samples (baseline at end of interval for now)

stimBLI = dsearchn(stimTimes',stimBL');
respBLI = dsearchn(respTimes',respBL');

allBetasStim(:,:,:) = allBetasStim(:,:,:) - nanmean(allBetasStim(:,:,stimBLI(1):stimBLI(2)),3);
allBetasResp(:,:,:) = allBetasResp(:,:,:) - nanmean(allBetasResp(:,:,respBLI(1):respBLI(2)),3);
allBetasTS(:,:,:) = allBetasTS(:,:,:) - nanmean(allBetasTS(:,:,(2):(2+tsBL)),3);

stimChannel = eeg_chaninds(EEG,{stimChannelString});
stim = squeeze(nanmean(allBetasStim(:,stimChannel,:),1));
stimCI = squeeze(nanstd(allBetasStim(:,stimChannel,:),[],1));
stimCI = tVal * stimCI ./ sqrt(length(participants));

respChannel = eeg_chaninds(EEG,{respChannelString});
resp = squeeze(nanmean(allBetasResp(:,respChannel,:),1));
respCI = squeeze(nanstd(allBetasResp(:,respChannel,:),[],1));
respCI = tVal * respCI ./ sqrt(length(participants));

tsChannel = eeg_chaninds(EEG,{tsChannelString});
ts = squeeze(nanmean(allBetasTS(:,tsChannel,:),1));
tsCI = squeeze(nanstd(allBetasTS(:,tsChannel,:),[],1));
tsCI = tVal * tsCI ./ sqrt(length(participants));

% Color map and scale (topographoes)
numContour = 5;
topoColour = cbrewer('div', 'RdBu', 12);
topoColour = flip(topoColour);
topoColour(topoColour > 1) = 1;
plotLineColour = topoColour(1,:);

plot_settings;

fig1 = makefigure(22,12);
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.1], [0.02 0.08], [0.06 0.04]);

subplot(2,3,1);
%plot(stimTimes,stim,'LineWidth',plotLineWidth,'Color',plotLineColour);
bl = boundedline(stimTimes,stim,stimCI,'cmap',plotColours(1,:));
bl.LineWidth = plotLineWidth;
title('Stimulus');

ax = gca;
ax.FontSize = plotFontSize;
ax.Box = 'off';
ax.YLabel.String = 'Voltage (\muV)';
ax.XLabel.String = 'Time (s)';
% % ax.YLim = yLimits;
ax.XLim = [stimTimes(1) stimTimes(end)];
text(ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))*0.05,ax.YLim(2),stimChannelString,'FontSize',plotFontSize);


hold on;
area(stimWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(stimWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');

% xlim([0.008 0.8]);

subplot(2,3,2);
% plot(respTimes,resp,'LineWidth',plotLineWidth,'Color',plotLineColour);
bl = boundedline(respTimes,resp,respCI,'cmap',plotColours(2,:));
bl.LineWidth = plotLineWidth;
switch whichTask
    case 'production'
        title('Response');
    case 'perception'
        title('Probe');
end
ax = gca;
ax.FontSize = plotFontSize;
ax.Box = 'off';
ax.YLabel.String = 'Voltage (\muV)';
ax.XLabel.String = 'Time (s)';
ax.XLim = [respTimes(1) respTimes(end)];
% ax.YLim = yLimits;
text(ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))*0.05,ax.YLim(2),respChannelString,'FontSize',plotFontSize);

hold on;
area(respWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(respWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');

% xlim([-1.2 -0.008]);
subplot(2,3,3);
% plot(ts,'LineWidth',plotLineWidth,'Color',plotLineColour);
bl = boundedline(1:length(ts),ts,tsCI,'cmap',plotColours(3,:));
bl.LineWidth = plotLineWidth;
title('Scaled');

ax = gca;
ax.FontSize = plotFontSize;
ax.Box = 'off';
ax.YLabel.String = 'Voltage (\muV)';
ax.XLabel.String = 'Interval proportion (%)';
ax.XLim = [1 scaleLength];
ax.XTick = [1 scaleLength/2 scaleLength];
ax.XAxis.TickLabels = {'1','50','100'};
% ax.YLim = [-300 100];
% ax.YLim = yLimits;
text(ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))*0.05,ax.YLim(2),tsChannelString,'FontSize',plotFontSize);


hold on;
area(scaledWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(scaledWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');


% Topographies
switch whichTask
    
    case {'production','perception'}
        mapLimits = [-8 8];
    case 'prediction'
        mapLimits = [-12 12];
end

stimWinPnts = dsearchn(stimTimes',stimWindow');
respWinPnts = dsearchn(respTimes',respWindow');
scaledWinPnts = scaledWindow(1):scaledWindow(2);

stimTopo = mean(mean(allBetasStim(:,:,stimWinPnts(1):stimWinPnts(2)),3),1);
respTopo = mean(mean(allBetasResp(:,:,respWinPnts(1):respWinPnts(2)),3),1);
scaledTopo = mean(mean(allBetasTS(:,:,scaledWinPnts(1):scaledWinPnts(2)),3),1);


subplot(2,3,4);
st = topoplot(stimTopo,EEG.chanlocs,'maplimits' ,mapLimits,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
st.Parent.XLim = [-0.6 0.6];
st.Parent.YLim = [-0.6 0.6];
c = colorbar;
c.Label.String = 'Voltage (\muV)';

subplot(2,3,5);
st = topoplot(respTopo,EEG.chanlocs,'maplimits' ,mapLimits,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
st.Parent.XLim = [-0.6 0.6];
st.Parent.YLim = [-0.6 0.6];

c = colorbar;
c.Label.String = 'Voltage (\muV)';

subplot(2,3,6);
st = topoplot(scaledTopo,EEG.chanlocs,'maplimits' ,mapLimits,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
st.Parent.XLim = [-0.6 0.6];
st.Parent.YLim = [-0.6 0.6];

c = colorbar;
c.Label.String = 'Voltage (\muV)';
% print(['../figures/fig_' whichTask '_rerps.tiff'],'-dtiff',tiffRes);

%% Sliding window
numWindows = 10;

subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.01], [0.02 0.08], [0.06 0.04]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 
fig = makefigure(36,12);
allTopoValues = [];
for w = 1:numWindows
    scaledWindow = round([(w-1)*scaleLength/10+1 (w)*scaleLength/10]);
    scaledWinPnts = scaledWindow(1):scaledWindow(2);
    
    % rERPs
    subplot(2,numWindows,w);
    bl = boundedline(1:length(ts),ts,tsCI,'cmap',plotLineColour);
    bl.LineWidth = plotLineWidth;
    titleString = [num2str(100*(scaledWindow(1)-1)/scaleLength) '-' num2str(100*scaledWindow(2)/scaleLength) '%'];
    title(titleString);
    
    ax = gca;
    ax.FontSize = plotFontSize;
    ax.Box = 'off';
%     ax.YLabel.String = 'Beta (a.u.)';
%     ax.XLabel.String = 'Proportion of interval (%)';
    % ax.XLim = [2 scaleLength-2];
%     ax.XTick = [1 scaleLength/2 scaleLength];
%     ax.XAxis.TickLabels = {'1','50','100'};
    % ax.YLim = [-300 100];
    % % ax.YLim = yLimits;
%     text(ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))*0.05,ax.YLim(2),tsChannelString,'FontSize',plotFontSize);
    
    
    hold on;
    area(scaledWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    area(scaledWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    
    % Topo
    scaledTopo = mean(mean(allBetasTS(:,:,scaledWinPnts(1):scaledWinPnts(2)),3),1);
    allTopoValues = [allTopoValues; scaledTopo];
end

switch whichTask
    
    case {'production','perception'}
        mapLimits = [-8 8];
    case 'prediction'
        mapLimits = [-12 12];
end

for w = 1:numWindows
    subplot(2,numWindows,numWindows+w);
    st = topoplot(allTopoValues(w,:),EEG.chanlocs,'maplimits' ,mapLimits,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
    st.Parent.XLim = [-0.6 0.6];
    st.Parent.YLim = [-0.6 0.6];
end


return;
%% Baseline correction/averaging

tVal = abs(tinv(0.025,length(participants)-1));

stimTimes = stimScalePnts/EEG.srate;
respTimes = respScalePnts/EEG.srate;

stimBLI = dsearchn(stimTimes',stimBL');
respBLI = dsearchn(respTimes',respBL');

stimChannelString = 'FCz';
respChannelString = 'FCz';

allBetasNoTSStim(:,:,:) = allBetasNoTSStim(:,:,:) - nanmean(allBetasNoTSStim(:,:,stimBLI(1):stimBLI(2)),3);
allBetasNoTSResp(:,:,:) = allBetasNoTSResp(:,:,:) - nanmean(allBetasNoTSResp(:,:,respBLI(1):respBLI(2)),3);

stimChannel = eeg_chaninds(EEG,{stimChannelString});
stim = squeeze(nanmean(allBetasNoTSStim(:,stimChannel,:),1));
stimCI = squeeze(nanstd(allBetasNoTSStim(:,stimChannel,:),[],1));
stimCI = tVal * stimCI ./ sqrt(length(participants));

respChannel = eeg_chaninds(EEG,{respChannelString});
resp = squeeze(nanmean(allBetasNoTSResp(:,respChannel,:),1));
respCI = squeeze(nanstd(allBetasNoTSResp(:,respChannel,:),[],1));
respCI = tVal * respCI ./ sqrt(length(participants));

% Analysis windows (topographies)
stimWindow = [0.084 0.124]; % N1
% respWindow = [-0.060 -0.008];
respWindow = [0 0.040];

% Color map and scale (topographoes)
mapMax = 0.5;
mapMin = -4.5;
numContour = 4;
%topoColour=cbrewer('seq', 'Blues', numContour);
topoColour=cbrewer('seq', 'Blues', numContour);
topoColour = flip(topoColour);

plotLineColour = topoColour(1,:);

plot_settings;

fig2 = makefigure(14,12);
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.1], [0.02 0.08], [0.06 0.04]);

subplot(2,2,1);
%plot(stimTimes,stim,'LineWidth',plotLineWidth,'Color',plotLineColour);
bl = boundedline(stimTimes,stim,stimCI,'cmap',plotLineColour);
bl.LineWidth = plotLineWidth;
title('Stimulus (N1)');

ax = gca;
ax.FontSize = plotFontSize;
ax.Box = 'off';
ax.YLabel.String = 'Voltage (\muV)';
ax.XLabel.String = 'Time (s)';
text(ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))*0.01,ax.YLim(2),stimChannelString,'FontSize',plotFontSize);

hold on;
area(stimWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(stimWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');

% xlim([0.008 0.8]);
subplot(2,2,2);
% plot(respTimes,resp,'LineWidth',plotLineWidth,'Color',plotLineColour);
bl = boundedline(respTimes,resp,respCI,'cmap',plotLineColour);
bl.LineWidth = plotLineWidth;
title('Response (BP)');

ax = gca;
ax.FontSize = plotFontSize;
ax.Box = 'off';
ax.YLabel.String = 'Voltage (\muV)';
ax.XLabel.String = 'Time (s)';
text(ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))*0.01,ax.YLim(2),respChannelString,'FontSize',plotFontSize);

hold on;
area(respWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(respWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');


% Topographies

stimWinPnts = dsearchn(stimTimes',stimWindow');
respWinPnts = dsearchn(respTimes',respWindow');
scaledWinPnts = scaledWindow(1):scaledWindow(2);

stimTopo = mean(mean(allBetasNoTSStim(:,:,stimWinPnts(1):stimWinPnts(2)),3),1);
respTopo = mean(mean(allBetasNoTSResp(:,:,respWinPnts(1):respWinPnts(2)),3),1);


subplot(2,2,3);
st = topoplot(stimTopo,EEG.chanlocs,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
st.Parent.XLim = [-0.6 0.6];
st.Parent.YLim = [-0.6 0.6];
c = colorbar;

subplot(2,2,4);
st = topoplot(respTopo,EEG.chanlocs,'numcontour' ,numContour,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',topoColour);
st.Parent.XLim = [-0.6 0.6];
st.Parent.YLim = [-0.6 0.6];

c = colorbar;

%print([../figures/whichTask '_rerp_with_no_scaled_component.tiff'],'-dtiff',tiffRes);
%% Error
grandAverageError = mean(allError,1);
% error1 = grandAverageError;
% error01 = grandAverageError;
% [h,iParticipant] = ttest(error1 - error01);

%% Residual
stimTimes = stimTimelimits(1):1/srate:stimTimelimits(2);
respTimes = respTimelimits(1):1/srate:respTimelimits(2);
clear subplot
whichChannel = 10;
whichChannel = 15;
meanStimR = squeeze(mean(allStimRs(:,:,whichChannel,:),1));
meanRespR = squeeze(mean(allRespRs(:,:,whichChannel,:),1));
figure;
subplot(1,2,1);
plot(stimTimes,meanStimR');
legend('Short','Medium','Long');
ylim([-1 1]);
title(['Stim ' num2str(mean(mean(meanStimR)))]);
subplot(1,2,2);
plot(respTimes,meanRespR');
legend('Short','Medium','Long');
ylim([-1 1]);
title(['Response ' num2str(mean(mean(meanRespR)))]);

whichChannel = 10;
whichChannel = 15;
meanStimRNoTS = squeeze(mean(allStimRsNoTS(:,:,whichChannel,:),1));
meanRespRNoTS = squeeze(mean(allRespRsNoTS(:,:,whichChannel,:),1));
figure;
subplot(1,2,1);
plot(stimTimes,meanStimRNoTS');
legend('Short','Medium','Long');
ylim([-1 1]);
title(['Stim (No TS) ' num2str(mean(mean(meanStimRNoTS)))]);
subplot(1,2,2);
plot(respTimes,meanRespRNoTS');
legend('Short','Medium','Long');
ylim([-1 1]);
title(['Response (No TS) ' num2str(mean(mean(meanRespRNoTS)))]);
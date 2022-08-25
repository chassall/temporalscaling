% Load preprocessed EEG and make ERPs
%
% Project: Temporal Scaling
% Other m-files required: EEGLAB, make_erp.m, find_artifacts.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% March 2020; Last revision: 12-Jul-2022

% Set data and results folders (change as needed)
dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task1-2_productionperception/data';
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

participants = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};
nParticipants = length(participants);
judgeTimes = [0.512 0.64 0.8 1 1.25; 1.056 1.32 1.65 2.0625 2.5781; 1.6 2 2.5 3.125 3.9063];
subconditions = {'very early','early','on time','late','very late'};
srate = 200;

%% Variables

% ERP windows
doShortBeepInterval = [-0.2 0.8];
doMediumBeepInterval = [-0.2 1.65];
doLongBeepInterval = [-0.2 2.5];
doShortResponseInterval = [-0.8 0.2] + [-0.2 0];
doMediumResponseInterval = [-1.652 0.2] + [-0.2 0];
doLongResponseInterval = [-2.5 0.2] + [-0.2 0];

judgeShortBeepInterval = [-0.2 0.8];
judgeMediumBeepInterval = [-0.2 1.652];
judgeLongBeepInterval = [-0.2 2.5];
judgeShortProbeInterval = [-0.8 0.2];
judgeMediumProbeInterval = [-1.652 0.2];
judgeLongProbeInterval = [-2.5 0.2];

% ERP baselines
doBeepBaseline = [-200 0];
doResponseBaseline = [-20 20];
judgeBeepBaseline = [-200 0];
judgeProbeBaseline = [-20 20];

% ERP triggers
allBeepTriggers = {12,22,32,42,52,62};
doShortBeepTriggers = {12};
doMediumBeepTriggers = {22};
doLongBeepTriggers = {32};
doShortResponseTriggers = {13};
doMediumResponseTriggers = {23};
doLongResponseTriggers = {33};
judgeShortBeepTriggers = {42};
judgeMediumBeepTriggers = {52};
judgeLongBeepTriggers = {62};
judgeShortProbeTriggers = {43};
judgeMediumProbeTriggers = {53};
judgeLongProbeTriggers = {63};

% ERP data
doShortBeepERP = nan(nParticipants, 32, round(200 * (doShortBeepInterval(2) - doShortBeepInterval(1))));
doMediumBeepERP = nan(nParticipants, 32, round(200 * (doMediumBeepInterval(2) - doMediumBeepInterval(1))));
doLongBeepERP = nan(nParticipants, 32, round(200 * (doLongBeepInterval(2) - doLongBeepInterval(1))));
doShortResponseERP = nan(nParticipants, 32, round(200 * (doShortResponseInterval(2) - doShortResponseInterval(1))));
doMediumResponseERP = nan(nParticipants, 32, round(200 * (doMediumResponseInterval(2) - doMediumResponseInterval(1))));
doLongResponseERP = nan(nParticipants, 32, round(200 * (doLongResponseInterval(2) - doLongResponseInterval(1))));
judgeShortBeepERP = nan(nParticipants, 32, round(200 * (judgeShortBeepInterval(2) - judgeShortBeepInterval(1))));
judgeMediumBeepERP = nan(nParticipants, 32, round(200 * (judgeMediumBeepInterval(2) - judgeMediumBeepInterval(1))));
judgeLongBeepERP = nan(nParticipants, 32, round(200 * (judgeLongBeepInterval(2) - judgeLongBeepInterval(1))));
judgeShortProbeERP = nan(nParticipants, 32, round(200 * (judgeShortProbeInterval(2) - judgeShortProbeInterval(1))));
judgeMediumProbeERP = nan(nParticipants, 32, round(200 * (judgeMediumProbeInterval(2) - judgeMediumProbeInterval(1))));
judgeLongProbeERP = nan(nParticipants, 32, round(200 * (judgeLongProbeInterval(2) - judgeLongProbeInterval(1))));

% Scaled ERP data - stretch or compress to match longest interval
doShortBeepERPS = nan(nParticipants, 32, round(200 * (doLongBeepInterval(2) - doLongBeepInterval(1))));
doMediumBeepERPS = nan(nParticipants, 32, round(200 * (doLongBeepInterval(2) - doLongBeepInterval(1))));
doLongBeepERPS = nan(nParticipants, 32, round(200 * (doLongBeepInterval(2) - doLongBeepInterval(1))));
doShortResponseERPS = nan(nParticipants, 32, round(200 * (doLongResponseInterval(2) - doLongResponseInterval(1))));
doMediumResponseERPS = nan(nParticipants, 32, round(200 * (doLongResponseInterval(2) - doLongResponseInterval(1))));
doLongResponseERPS = nan(nParticipants, 32, round(200 * (doLongResponseInterval(2) - doLongResponseInterval(1))));
judgeShortBeepERPS = nan(nParticipants, 32, round(200 * (judgeLongBeepInterval(2) - judgeLongBeepInterval(1))));
judgeMediumBeepERPS = nan(nParticipants, 32, round(200 * (judgeLongBeepInterval(2) - judgeLongBeepInterval(1))));
judgeLongBeepERPS = nan(nParticipants, 32, round(200 * (judgeLongBeepInterval(2) - judgeLongBeepInterval(1))));
judgeShortProbeERPS = nan(nParticipants, 32, round(200 * (judgeLongProbeInterval(2) - judgeLongProbeInterval(1))));
judgeMediumProbeERPS = nan(nParticipants, 32, round(200 * (judgeLongProbeInterval(2) - judgeLongProbeInterval(1))));
judgeLongProbeERPS = nan(nParticipants, 32, round(200 * (judgeLongProbeInterval(2) - judgeLongProbeInterval(1))));

% Artifact proportions for each ERP of interest (see header for order)
allArtifacts = nan(nParticipants, 12, 32);

% Loop through participants
for iParticipant = 1:nParticipants
    
    % Load preprocessed data
    subString = ['sub-' participants{iParticipant}];
    preprocessedFolder = fullfile(dataFolder,'derivatives','eegprep',subString);
    preprocessedFile = [subString '_task-temporalscaling_eegprep'];
    load(fullfile(preprocessedFolder,preprocessedFile), 'EEG');
    
    % Remove ocular components using the results of ICLabel
    eyeLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Eye'));
    [~,I] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
    whichOnes = find(I == eyeLabel); % Max out of all possible labels
    EEG = pop_subcomp(EEG,whichOnes,0);
    
    % Test ERP (check for back channelss
    [~, ~,~,~,theseArtifactsByChannel] = make_erp(EEG,allBeepTriggers,doShortBeepInterval,doBeepBaseline);
    artifactProps = mean(theseArtifactsByChannel,2);
    if any(artifactProps > 0.2)
        error('bad channel found');
        break;
    end

    % Make ERPs
    artifactSettings.maxMin = 150;
    artifactSettings.level = 150;
    artifactSettings.step = 50;
    artifactSettings.lowest = 0.1;

    [doShortBeepERP(iParticipant,:,:), allArtifacts(iParticipant,1,:),doShortBeepTimes,~,~] = make_erp(EEG,doShortBeepTriggers,doShortBeepInterval,doBeepBaseline,artifactSettings);
    [doMediumBeepERP(iParticipant,:,:), allArtifacts(iParticipant,2,:),doMediumBeepTimes] = make_erp(EEG,doMediumBeepTriggers,doMediumBeepInterval,doBeepBaseline);
    [doLongBeepERP(iParticipant,:,:), allArtifacts(iParticipant,3,:),doLongBeepTimes] = make_erp(EEG,doLongBeepTriggers,doLongBeepInterval,doBeepBaseline);  
    [doShortResponseERP(iParticipant,:,:), allArtifacts(iParticipant,4,:),doShortResponseTimes] = make_erp(EEG,doShortResponseTriggers,doShortResponseInterval,doResponseBaseline);
    [doMediumResponseERP(iParticipant,:,:), allArtifacts(iParticipant,5,:),doMediumResponseTimes] = make_erp(EEG,doMediumResponseTriggers,doMediumResponseInterval,doResponseBaseline);
    [doLongResponseERP(iParticipant,:,:), allArtifacts(iParticipant,6,:),doLongResponseTimes] = make_erp(EEG,doLongResponseTriggers,doLongResponseInterval,doResponseBaseline); 
    [judgeShortBeepERP(iParticipant,:,:), allArtifacts(iParticipant,7,:),judgeShortBeepTimes] = make_erp(EEG,judgeShortBeepTriggers,judgeShortBeepInterval,judgeBeepBaseline);
    [judgeMediumBeepERP(iParticipant,:,:), allArtifacts(iParticipant,8,:),judgeMediumBeepTimes] = make_erp(EEG,judgeMediumBeepTriggers,judgeMediumBeepInterval,judgeBeepBaseline);
    [judgeLongBeepERP(iParticipant,:,:), allArtifacts(iParticipant,9,:),judgeLongBeepTimes] = make_erp(EEG,judgeLongBeepTriggers,judgeLongBeepInterval,judgeBeepBaseline);  
    [judgeShortProbeERP(iParticipant,:,:), allArtifacts(iParticipant,10,:),judgeShortProbeTimes] = make_erp(EEG,judgeShortProbeTriggers,judgeShortProbeInterval,judgeProbeBaseline);
    [judgeMediumProbeERP(iParticipant,:,:), allArtifacts(iParticipant,11,:),judgeMediumProbeTimes] = make_erp(EEG,judgeMediumProbeTriggers,judgeMediumProbeInterval,judgeProbeBaseline);
    [judgeLongProbeERP(iParticipant,:,:), allArtifacts(iParticipant,12,:),judgeLongProbeTimes] = make_erp(EEG,judgeLongProbeTriggers,judgeLongProbeInterval,judgeProbeBaseline); 


end

% Grand averages
doShortBeepGrandAve = squeeze(nanmean(doShortBeepERP,1));
doMediumBeepGrandAve = squeeze(nanmean(doMediumBeepERP,1));
doLongBeepGrandAve = squeeze(nanmean(doLongBeepERP,1));
doShortResponseGrandAve = squeeze(nanmean(doShortResponseERP,1));
doMediumResponseGrandAve = squeeze(nanmean(doMediumResponseERP,1));
doLongResponseGrandAve = squeeze(nanmean(doLongResponseERP,1));
judgeShortBeepGrandAve = squeeze(nanmean(judgeShortBeepERP,1));
judgeMediumBeepGrandAve = squeeze(nanmean(judgeMediumBeepERP,1));
judgeLongBeepGrandAve = squeeze(nanmean(judgeLongBeepERP,1));
judgeShortProbeGrandAve = squeeze(nanmean(judgeShortProbeERP,1));
judgeMediumProbeGrandAve = squeeze(nanmean(judgeMediumProbeERP,1));
judgeLongProbeGrandAve = squeeze(nanmean(judgeLongProbeERP,1));

% Save everything
save(fullfile(resultsFolder,'erp_results_production_and_perception_tasks.mat'));
return;

%% Unused in manuscript (only for if a participant requests to see their brain activity)

channelS = 'FCz';
channelI = eeg_chaninds(EEG,channelS);

makefigure(33.87,19.05);

axs = {};

axs{1} = subplot(2,2,1);
plot(doShortBeepTimes,doShortBeepGrandAve(channelI,:));
hold on;
plot(doMediumBeepTimes,doMediumBeepGrandAve(channelI,:));
plot(doLongBeepTimes,doLongBeepGrandAve(channelI,:));
title('''Do'' Task - Beep');

axs{2} = subplot(2,2,2);
plot(doShortResponseTimes,doShortResponseGrandAve(channelI,:)); hold on;
plot(doMediumResponseTimes,doMediumResponseGrandAve(channelI,:));
plot(doLongResponseTimes,doLongResponseGrandAve(channelI,:));
title('''Do'' Task - Response');

axs{3} = subplot(2,2,3);
plot(judgeShortBeepTimes,judgeShortBeepGrandAve(channelI,:));
hold on;
plot(judgeMediumBeepTimes,judgeMediumBeepGrandAve(channelI,:));
plot(judgeLongBeepTimes,judgeLongBeepGrandAve(channelI,:));
title('''Judge'' Task - First Beep');

axs{4} = subplot(2,2,4);
plot(judgeShortProbeTimes,judgeShortProbeGrandAve(channelI,:));
hold on;
plot(judgeMediumProbeTimes,judgeMediumProbeGrandAve(channelI,:));
plot(judgeLongProbeTimes,judgeLongProbeGrandAve(channelI,:));
title('''Judge'' Task - Second Beep');

for i = 1:4
    axs{i}.XLabel.String = 'Time (ms)';
    axs{i}.YLabel.String = 'Voltage (\muV)';
    axs{i}.Box = 'off';
    axs{i}.FontSize = 12;
    legend(axs{i},{'short','medium','long'},'Box','off','Location','SouthWest');

    for j = 1:length(axs{i}.Children)
        axs{i}.Children(j).LineWidth = 1.5;
    end
end

print('scratch_erps.jpg','-djpeg','-r300');
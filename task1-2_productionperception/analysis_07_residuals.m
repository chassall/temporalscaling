% Make ERPs based on the residual EEG
%
% Project: Temporal Scaling
% Other m-files required: EEGLAB, make_erp.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% March 2020; Last revision: 12-Jul-2022

close all; clear all; clc; rng(41); % for reproducibility

% Set results folders (change as needed)
dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task1-2_productionperception/data';
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

participants = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};
nParticipants = length(participants);
whichTask = 'production';
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

% Load ERPs
load(fullfile(resultsFolder,'erp_results_production_and_perception_tasks.mat'));

targetTimes = [0.8 1.65 2.5];
doShortBeepInterval = [-0.2 targetTimes(1) + 0.2];
doMediumBeepInterval = [-0.2 targetTimes(2) + 0.2];
doLongBeepInterval = [-0.2 targetTimes(3) + 0.2];

for iParticipant = 1:nParticipants

    subjString = ['sub-' participants{iParticipant}];
    disp(subjString);
    
    % Load preprocessed data
    preprocessedFolder = fullfile(dataFolder,'derivatives','eegprep',subjString);
    preprocessedFile = [subjString '_task-temporalscaling_eegprep'];
    load(fullfile(preprocessedFolder,preprocessedFile), 'EEG');
    
    % Remove ocular components using the results of ICLabel
    eyeLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Eye'));
    [~,I] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
    whichOnes = find(I == eyeLabel); % Max out of all possible labels
    EEG = pop_subcomp(EEG,whichOnes,0);
    
    % Load the residuals, which will be in EEGLAB format (called rEEG)
    residualsFolder = fullfile(dataFolder,'derivatives','eegresidual',subjString);
    residualsFile = [subjString '_task-' whichTask '_eegresidualbox'];
    load(fullfile(residualsFolder,residualsFile),'rEEGWithTS','rEEGWithFixed','tsBeta','toExclude');
    
    % First, ERPs only
    [doShortERP(iParticipant,:,:)] = make_erp(EEG,doShortBeepTriggers,doShortBeepInterval,doBeepBaseline);
    [doMediumERP(iParticipant,:,:)] = make_erp(EEG,doMediumBeepTriggers,doMediumBeepInterval,doBeepBaseline);
    [doLongERP(iParticipant,:,:)] = make_erp(EEG,doLongBeepTriggers,doLongBeepInterval,doBeepBaseline);  
    
    % First, make ERPs with fixed only
    [doShortERPF(iParticipant,:,:)] = make_erp(rEEGWithFixed,doShortBeepTriggers,doShortBeepInterval,doBeepBaseline);
    [doMediumERPF(iParticipant,:,:)] = make_erp(rEEGWithFixed,doMediumBeepTriggers,doMediumBeepInterval,doBeepBaseline);
    [doLongERPF(iParticipant,:,:)] = make_erp(rEEGWithFixed,doLongBeepTriggers,doLongBeepInterval,doBeepBaseline);  

    % Now scaled only
    [doShortERPS(iParticipant,:,:)] = make_erp(rEEGWithTS,doShortBeepTriggers,doShortBeepInterval,doBeepBaseline);
    [doMediumERPS(iParticipant,:,:)] = make_erp(rEEGWithTS,doMediumBeepTriggers,doMediumBeepInterval,doBeepBaseline);
    [doLongERPS(iParticipant,:,:)] = make_erp(rEEGWithTS,doLongBeepTriggers,doLongBeepInterval,doBeepBaseline);  
   

end

save(fullfile(resultsFolder,'residual_results_production_task.mat'));

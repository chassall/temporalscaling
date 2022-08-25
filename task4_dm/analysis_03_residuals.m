% Residual analysis in line with rest of temporal scaling paper

% Load preprocessed EEG and make ERPs
%
% Project: Temporal Scaling
% Other m-files required: EEGLAB, make_erp.m, find_artifacts.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% March 2020; Last revision: 02-Dec-2020

rawFolder = 'E:\ccn_lab_raw\ds002734-download\';
participants = {'02','03','04','05','06','07','08','10','11','12','13','14','15','16','18','19','20','21'};
nParticipants = length(participants);
judgeTimes = [0.512 0.64 0.8 1 1.25; 1.056 1.32 1.65 2.0625 2.5781; 1.6 2 2.5 3.125 3.9063];
subconditions = {'fast','slow'};
srate = 200;
nElectrodes = 62;
doLaplacian = 0;
%% Variables

% ERP windows
stimInterval = [-0.2,1.5];

% ERP baselines
stimBaseline = [-200,0];

% ERP triggers
fastStimTriggers = {3};
slowStimTriggers = {4};
fastRespTriggers = {5};
slowRespTriggers = {6};

% ERP data
fastStimERP = nan(nParticipants,nElectrodes, round(srate * (stimInterval(2) - stimInterval(1))));
slowStimERP = nan(nParticipants,nElectrodes, round(srate * (stimInterval(2) - stimInterval(1))));
fastStimERPF = nan(nParticipants,nElectrodes, round(srate * (stimInterval(2) - stimInterval(1))));
slowStimERPF = nan(nParticipants,nElectrodes, round(srate * (stimInterval(2) - stimInterval(1))));
fastStimERPS = nan(nParticipants,nElectrodes, round(srate * (stimInterval(2) - stimInterval(1))));
slowStimERPS = nan(nParticipants,nElectrodes, round(srate * (stimInterval(2) - stimInterval(1))));

% Artifact proportions for each ERP of interest (see header for order)
allArtifacts = nan(nParticipants, 4);

meanRTs = [];

% Loop through participants
for iParticipant = 1:nParticipants
    
    % Load preprocessed data
    load(['sub-' participants{iParticipant}]);
    EEG = pop_reref(EEG,{'TP7','TP8'});
    
    % Load the residuals, which will be in EEGLAB format (called rEEG)
    residualsFolder = ['./data/derivatives/eegresidual/sub-' participants{iParticipant}];
    residualsFile = ['sub-' participants{iParticipant} '_task-dm_eegresidual_reref'];
    load(fullfile(residualsFolder,residualsFile),'rEEGWithTS','rEEGWithFixed','beta','toExclude');
    
    % Do a media split on RTs
    % Get RTs from EEG
    eventTimes = [EEG.event.latency];
    eventRTs = eventTimes(2:2:end) - eventTimes(1:2:end-1);
    meanRTs(iParticipant) = mean(eventRTs);
    isSlow = eventRTs > median(eventRTs);
    
    % Set event names/conditions
    trialCount = 1;
    for i = 1:length(EEG.event)
       if EEG.event(i).type == 1
           if isSlow(trialCount)
               EEG.event(i).type = 4;
           else
               EEG.event(i).type = 3;
           end
       elseif EEG.event(i).type == 2
           if isSlow(trialCount)
               EEG.event(i).type = 6;
           else
               EEG.event(i).type = 5;
           end
           trialCount = trialCount + 1;
       end
    end
    
    rEEGWithFixed.event = EEG.event;
    rEEGWithTS.event = EEG.event;
    
    % Make ERPs
    [fastStimERP(iParticipant,:,:), allArtifacts(iParticipant,1)] = make_erp(EEG,fastStimTriggers,stimInterval,stimBaseline);
    [slowStimERP(iParticipant,:,:), allArtifacts(iParticipant,2)] = make_erp(EEG,slowStimTriggers,stimInterval,stimBaseline);
    [fastStimERPF(iParticipant,:,:), allArtifactsF(iParticipant,1)] = make_erp(rEEGWithFixed,fastStimTriggers,stimInterval,stimBaseline);
    [slowStimERPF(iParticipant,:,:), allArtifactsF(iParticipant,2)] = make_erp(rEEGWithFixed,slowStimTriggers,stimInterval,stimBaseline);
    [fastStimERPS(iParticipant,:,:), allArtifactsS(iParticipant,1)] = make_erp(rEEGWithTS,fastStimTriggers,stimInterval,stimBaseline);
    [slowStimERPS(iParticipant,:,:), allArtifactsS(iParticipant,2)] = make_erp(rEEGWithTS,slowStimTriggers,stimInterval,stimBaseline);
end

%%
% Grand averages
fastStimGrandAve = squeeze(nanmean(fastStimERP));
slowStimGrandAve = squeeze(nanmean(slowStimERP));
fastStimGrandAveF = squeeze(nanmean(fastStimERPF));
slowStimGrandAveF = squeeze(nanmean(slowStimERPF));
fastStimGrandAveS = squeeze(nanmean(fastStimERPS));
slowStimGrandAveS = squeeze(nanmean(slowStimERPS));

%%
save('residuals.mat');

%%
load('residuals.mat');
iElectrode = 31;
axs = [];
axs(1) = subplot(1,3,1); plot(fastStimGrandAve(iElectrode,:)); hold on;  plot(slowStimGrandAve(iElectrode,:)); title('all'); legend('fast','slow');
axs(2) = subplot(1,3,2); plot(fastStimGrandAveF(iElectrode,:)); hold on;  plot(slowStimGrandAveF(iElectrode,:)); title('fixed');
axs(3) = subplot(1,3,3); plot(fastStimGrandAveS(iElectrode,:)); hold on;  plot(slowStimGrandAveS(iElectrode,:)); title('scaled');

linkaxes(axs);
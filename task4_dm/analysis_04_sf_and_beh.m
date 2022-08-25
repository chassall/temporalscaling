% Scale factor analysis for the DM task
%
% Project: Temporal Scaling

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% August 2022

close all; clear all; clc;
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

% Compute scale factors

participants = {'02','03','04','05','06','07','08','10','11','12','13','14','15','16','18','19','20','21'};
nParticipants = length(participants);

meanRT = 0.765;
dmTimes = [meanRT]; % Mean decision time
dmPnts = dmTimes * 200; % Multiply by sample rate

% Rhythmic and repeated
whichConditions = {{131,231},{142,242},{133,233},{144,244}}; % s cue, l cue, s target, l target

allRPre = [];
allRScaled = [];
allRFixed = [];
allR2Pre = [];
allR2Scaled = [];
allR2Fixed = [];

allRTs = [];

% To compute a grand average
combinedSignalsPre = [];
combinedSignalsFixed= [];
combinedSignalsScaled = [];


iElectrode = 35; % C1, where the scaled signal was maximal

bl = [];

% Loop through participants
for iParticipant = 1:nParticipants
    
    % Load preprocessed data
    preprocessedFolder = './';
    preprocessedFile = ['sub-' participants{iParticipant}];
    load(fullfile(preprocessedFolder,preprocessedFile), 'EEG');

    % Load residual data
    preprocessedFolder = ['./data/derivatives/eegresidual/sub-' participants{iParticipant}];
    preprocessedFile = ['sub-' participants{iParticipant} '_task-dm_eegresidual_reref'];
    load(fullfile(preprocessedFolder,preprocessedFile),'rEEGWithFixed','rEEGWithTS','toExclude');
    
    
    % Loop through events, gather epochs spanning cue to response, and
    % stretch as needed
    preEpochsS = [];
    preEpochsL = [];
    fixedEpochsS = [];
    fixedEpochsL = [];
    scaledEpochsS = [];
    scaledEpochsL = [];
    for i = 1:length(EEG.event)-1
        
        % Is this a cue?
        if EEG.event(i).type == 1
            
            preEpoch = EEG.data(:,EEG.event(i).latency:EEG.event(i+1).latency);
            fixedEpoch = rEEGWithFixed.data(:,EEG.event(i).latency:EEG.event(i+1).latency);
            scaledEpoch = rEEGWithTS.data(:,EEG.event(i).latency:EEG.event(i+1).latency);
            
                
                thisRT = size(preEpoch,2) * (1/EEG.srate);
                allRTs = [allRTs thisRT];
                
                % Short trial
                if thisRT <= meanRT
                    preEpochsS(:,:,end+1) = imresize(preEpoch,[64,dmPnts(1)]);
                    fixedEpochsS(:,:,end+1) = imresize(fixedEpoch,[64,dmPnts(1)]);
                    scaledEpochsS(:,:,end+1) = imresize(scaledEpoch,[64,dmPnts(1)]);
                % Long trial
                else
                    preEpochsL(:,:,end+1) = imresize(preEpoch,[64,dmPnts(1)]);
                    fixedEpochsL(:,:,end+1) = imresize(fixedEpoch,[64,dmPnts(1)]);
                    scaledEpochsL(:,:,end+1) = imresize(scaledEpoch,[64,dmPnts(1)]);
                end
            
        end
        
    end
    

    maxMin = 150;
    level = 150;
    step = 40;
    lowest = 0.1;

    dEEG = EEG;
    dEEG.data = preEpochsS(:,:,:); dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanPreS = squeeze(mean(dEEG.data(iElectrode,:,~isArtifact),3));
    dEEG.data = preEpochsL(:,:,:); dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanPreL = squeeze(mean(dEEG.data(iElectrode,:,~isArtifact),3));
    allMeanPre = [meanPreS' meanPreL'];

    dEEG.data = fixedEpochsS(:,:,:); dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanFixedS = squeeze(mean(dEEG.data(iElectrode,:,~isArtifact),3));
    dEEG.data = fixedEpochsL(:,:,:); dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanFixedL = squeeze(mean(dEEG.data(iElectrode,:,~isArtifact),3));
    allMeanFixed = [meanFixedS' meanFixedL'];

    dEEG.data = scaledEpochsS(:,:,:); dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanScaledS = squeeze(mean(dEEG.data(iElectrode,:,~isArtifact),3));
    dEEG.data = scaledEpochsL(:,:,:); dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanScaledL = squeeze(mean(dEEG.data(iElectrode,:,~isArtifact),3));
    allMeanScaled = [meanScaledS' meanScaledL'];

    allRPre(iParticipant,:,:) = corr(allMeanPre);
    allRFixed(iParticipant,:,:) = corr(allMeanFixed);
    allRScaled(iParticipant,:,:) = corr(allMeanScaled);

    % Compute r-squared
    % Try to predict 'i' using 'j'
    for i = 1:2
        for j = 1:2
            allR2Pre(iParticipant,i,j) = myrsquare(allMeanPre(:,i),allMeanPre(:,j));
            allR2Fixed(iParticipant,i,j) = myrsquare(allMeanFixed(:,i),allMeanFixed(:,j));
            allR2Scaled(iParticipant,i,j) = myrsquare(allMeanScaled(:,i),allMeanScaled(:,j));
        end
    end


    combinedSignalsPre(iParticipant,:,:) = allMeanPre;
    combinedSignalsFixed(iParticipant,:,:) = allMeanFixed;
    combinedSignalsScaled(iParticipant,:,:) = allMeanScaled;

end

save(fullfile(resultsFolder,'scalefactor_dm.mat'));
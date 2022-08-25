% Compute scale factors in the Prediction task
%
% Other m-files required: eeglab

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% 25-Aug-2022

close all; clear all; clc;
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

ps = [1:12 14 16:21]; % 13, 15 missing on datadryad
participants = num2str(ps,'%0.2d,');
participants(end) = [];
participants = strsplit(participants,',');
nParticipants = length(participants);

predictionTimes = [0.7 1.3];
predictionPnts = predictionTimes * 200; % Multiply by sample rate

% Rhythmic and repeated
whichConditions = {{131,231},{142,242},{133,233},{144,244}}; % s cue, l cue, s target, l target

allRPre = [];
allRScaled = [];
allRFixed = [];
allR2Pre = [];
allR2Scaled = [];
allR2Fixed = [];

allRTSRhythmic = [];
allRTLRhythmic = [];
allRTSRepeated = [];
allRTLRepeated = [];

% To compute a grand average
combinedSignalsPre = [];
combinedSignalsFixed= [];
combinedSignalsScaled = [];

bl = [];

iElectrode = 47; % FCz, where the scaled signal was maximal

% Loop through participants
for iParticipant = 1:nParticipants
    
    % Load preprocessed data
    preprocessedFolder = ['./data/derivatives/eegprep/sub-' participants{iParticipant}];
    preprocessedFile = ['sub-' participants{iParticipant} '_task-entrainrhythm1b_eegprep'];
    load(fullfile(preprocessedFolder,preprocessedFile), 'EEG');

    % Load residual data
    preprocessedFolder = ['./data/derivatives/eegresidual/sub-' participants{iParticipant}];
    preprocessedFile = ['sub-' participants{iParticipant} '_task-prediction_eegresidualbox'];
    load(fullfile(preprocessedFolder,preprocessedFile),'rEEGWithFixed','rEEGWithTS','toExclude');
    
    
    tempRTSRhythmic = [];
    tempRTLRhythmic = [];
    tempRTSRepeated = [];
    tempRTLRepeated = [];
    
    % Loop through events, gather epochs spanning cue to response, and
    % stretch as needed
    preEpochsS = [];
    preEpochsL = [];
    fixedEpochsS = [];
    fixedEpochsL = [];
    scaledEpochsS = [];
    scaledEpochsL = [];
    for i = 1:length(EEG.event)-2
        
        % Is this a cue?
        if ismember(EEG.event(i).type,[131,231,142,242])
            
            preEpoch = EEG.data(:,EEG.event(i).latency:EEG.event(i+1).latency);
            fixedEpoch = rEEGWithFixed.data(:,EEG.event(i).latency:EEG.event(i+1).latency);
            scaledEpoch = rEEGWithTS.data(:,EEG.event(i).latency:EEG.event(i+1).latency);
            
            thisRT = (EEG.event(i+2).latency - EEG.event(i+1).latency) * (1/EEG.srate);
            
            switch EEG.event(i).type
                case 131
                    tempRTSRhythmic = [tempRTSRhythmic thisRT];
                case 231
                    tempRTLRhythmic = [tempRTLRhythmic thisRT];
                case 142
                    tempRTSRepeated = [tempRTSRepeated thisRT];
                case 242
                    tempRTLRepeated = [tempRTLRepeated thisRT];
            end
            
            
            switch EEG.event(i).type
                
                % Short cue
                case whichConditions{1}
                    preEpochsS(:,:,end+1) = imresize(preEpoch,[64,predictionPnts(2)]);
                    fixedEpochsS(:,:,end+1) = imresize(fixedEpoch,[64,predictionPnts(2)]);
                    scaledEpochsS(:,:,end+1) = imresize(scaledEpoch,[64,predictionPnts(2)]);
                % Long cue
                case whichConditions{2}
                    preEpochsL(:,:,end+1) = imresize(preEpoch,[64,predictionPnts(2)]);
                    fixedEpochsL(:,:,end+1) = imresize(fixedEpoch,[64,predictionPnts(2)]);
                    scaledEpochsL(:,:,end+1) = imresize(scaledEpoch,[64,predictionPnts(2)]);
            end
            
        end
        
    end
    

    allRTSRhythmic(iParticipant) = mean(tempRTSRhythmic);
    allRTLRhythmic(iParticipant) = mean(tempRTLRhythmic);
    allRTSRepeated(iParticipant) = mean(tempRTSRepeated);
    allRTLRepeated(iParticipant) = mean(tempRTLRepeated);
    
    maxMin = 150;
    level = 150;
    step = 40;
    lowest = 0.1;

    dEEG = EEG;
    dEEG.data = preEpochsS(:,:,:); 
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanPreS = squeeze(mean(dEEG.data(iElectrode,:,~isArtifact),3));
    dEEG.data = preEpochsL(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanPreL = squeeze(mean(dEEG.data(iElectrode,:,~isArtifact),3));
    allMeanPre = [meanPreS' meanPreL'];

    dEEG.data = fixedEpochsS(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanFixedS = squeeze(mean(dEEG.data(iElectrode,:,~isArtifact),3));
    dEEG.data = fixedEpochsL(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanFixedL = squeeze(mean(dEEG.data(iElectrode,:,~isArtifact),3));
    allMeanFixed = [meanFixedS' meanFixedL'];

    dEEG.data = scaledEpochsS(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanScaledS = squeeze(mean(dEEG.data(iElectrode,:,~isArtifact),3));
    dEEG.data = scaledEpochsL(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
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

save(fullfile(resultsFolder,'scalefactor_prediction.mat'));
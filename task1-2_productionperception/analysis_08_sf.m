% Compute scale factors
%
% Project: Temporal Scaling
% Other m-files required: EEGLAB

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% March 2020; Last revision: 12-Jul-2022

close all; clear all; clc;

% Set results folders (change as needed)
dataFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\task1-2_productionperception\data';
resultsFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\results';

whichTask = 'perception'; % 'production' or 'perception'

% Compute scale factors
participants = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};
nParticipants = length(participants);

productionTimes = [0.8 1.65 2.5];
productionPnts = productionTimes * 200; % Multiply by sample rate

allRPre = [];
allRScaled = [];
allRFixed = [];

allR2Pre = [];
allR2Scaled = [];
allR2Fixed = [];

bl = [];

% Loop through participants
for iParticipant = 1:nParticipants
    
    subjString = ['sub-' participants{iParticipant}];

    % Load preprocessed data
    preprocessedFolder = fullfile(dataFolder,'derivatives','eegprep',subjString);
    preprocessedFile = [subjString '_task-temporalscaling_eegprep'];
    load(fullfile(preprocessedFolder,preprocessedFile), 'EEG');

    % Load residual data
    preprocessedFolder = fullfile(dataFolder,'derivatives','eegresidual',subjString);
    preprocessedFile = [subjString '_task-' whichTask '_eegresidualbox'];
    load(fullfile(preprocessedFolder,preprocessedFile),'rEEGWithFixed','rEEGWithTS','toExclude');
    
    % EEG = rEEGWithFixed; % Pick one of the two residual data sets

    % Add labels to event variable
    for iEvent = 1:length(EEG.event)
        switch EEG.event(iEvent).type
            case 2
                EEG.event(iEvent).task = 'metronome';
                EEG.event(iEvent).type = 'beep';
                EEG.event(iEvent).condition = 'metronome';
            case 12
                EEG.event(iEvent).task = 'production';
                EEG.event(iEvent).type = 'beep';
                EEG.event(iEvent).condition = 'short';
            case 13
                EEG.event(iEvent).task = 'production';
                EEG.event(iEvent).type = 'response';
                EEG.event(iEvent).condition = 'short';
            case 22
                EEG.event(iEvent).task = 'production';
                EEG.event(iEvent).type = 'beep';
                EEG.event(iEvent).condition = 'medium';
            case 23
                EEG.event(iEvent).task = 'production';
                EEG.event(iEvent).type = 'response';
                EEG.event(iEvent).condition = 'medium';
            case 32
                EEG.event(iEvent).task = 'production';
                EEG.event(iEvent).type = 'beep';
                EEG.event(iEvent).condition = 'long';
            case 33
                EEG.event(iEvent).task = 'production';
                EEG.event(iEvent).type = 'response';
                EEG.event(iEvent).condition = 'long';
            case 42
                EEG.event(iEvent).task = 'perception';
                EEG.event(iEvent).type = 'beep';
                EEG.event(iEvent).condition = 'short';
            case 43
                EEG.event(iEvent).task = 'perception';
                EEG.event(iEvent).type = 'response';
                EEG.event(iEvent).condition = 'short';
            case 52
                EEG.event(iEvent).task = 'perception';
                EEG.event(iEvent).type = 'beep';
                EEG.event(iEvent).condition = 'medium';
            case 53
                EEG.event(iEvent).task = 'perception';
                EEG.event(iEvent).type = 'response';
                EEG.event(iEvent).condition = 'medium';
            case 62
                EEG.event(iEvent).task = 'perception';
                EEG.event(iEvent).type = 'beep';
                EEG.event(iEvent).condition = 'long';
            case 63
                EEG.event(iEvent).task = 'perception';
                EEG.event(iEvent).type = 'response';
                EEG.event(iEvent).condition = 'long';
            otherwise % No change
                EEG.event(iEvent).type = num2str(EEG.event(iEvent).type);
        end
    end

    % Loop through events, gather epochs spanning cue to response, and
    % stretch as needed
    preEpochsS = [];
    preEpochsM = [];
    preEpochsL = [];
    fixedEpochsS = [];
    fixedEpochsM = [];
    fixedEpochsL = [];
    scaledEpochsS = [];
    scaledEpochsM = [];
    scaledEpochsL = [];
    for i = 1:length(rEEGWithTS.event)
        if (strcmp(whichTask,'production') && strcmp(EEG.event(i).task,'production') && strcmp(EEG.event(i).type,'beep')) || ...
                (strcmp(whichTask,'perception') && strcmp(EEG.event(i).task,'perception') && strcmp(EEG.event(i).type,'beep'))
            
            preEpoch = EEG.data(:,EEG.event(i).latency:EEG.event(i+1).latency);
            fixedEpoch = rEEGWithFixed.data(:,EEG.event(i).latency:EEG.event(i+1).latency);
            scaledEpoch = rEEGWithTS.data(:,EEG.event(i).latency:EEG.event(i+1).latency);

            switch EEG.event(i).condition
                case 'short'
                    preEpochsS(:,:,end+1) = imresize(preEpoch,[32,productionPnts(3)]);
                    fixedEpochsS(:,:,end+1) = imresize(fixedEpoch,[32,productionPnts(3)]);
                    scaledEpochsS(:,:,end+1) = imresize(scaledEpoch,[32,productionPnts(3)]);

                    % sEpochs(:,:,end+1) = imresize(thisEpoch,[32,productionPnts(3)]);
                case 'medium'
                    preEpochsM(:,:,end+1) = imresize(preEpoch,[32,productionPnts(3)]);
                    fixedEpochsM(:,:,end+1) = imresize(fixedEpoch,[32,productionPnts(3)]);
                    scaledEpochsM(:,:,end+1) = imresize(scaledEpoch,[32,productionPnts(3)]);

                    % mEpochs(:,:,end+1) = imresize(thisEpoch,[32,productionPnts(3)]);
                case 'long'
                    preEpochsL(:,:,end+1) = imresize(preEpoch,[32,productionPnts(3)]);
                    fixedEpochsL(:,:,end+1) = imresize(fixedEpoch,[32,productionPnts(3)]);
                    scaledEpochsL(:,:,end+1) = imresize(scaledEpoch,[32,productionPnts(3)]);

                    % lEpochs(:,:,end+1) = imresize(thisEpoch,[32,productionPnts(3)]);
            end
        end
    end
    

    maxMin = 150;
    level = 150;
    step = 40;
    lowest = 0.1;

    dEEG = EEG;
    dEEG.data = preEpochsS(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanPreS = squeeze(mean(dEEG.data(15,:,~isArtifact),3));
   
    dEEG.data = preEpochsM(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanPreM = squeeze(mean(dEEG.data(15,:,~isArtifact),3));
    
    dEEG.data = preEpochsL(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanPreL = squeeze(mean(dEEG.data(15,:,~isArtifact),3));
    
    allMeanPre = [meanPreS' meanPreM' meanPreL'];

    dEEG.data = fixedEpochsS(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanFixedS = squeeze(mean(dEEG.data(15,:,~isArtifact),3));
    
    dEEG.data = fixedEpochsM(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanFixedM = squeeze(mean(dEEG.data(15,:,~isArtifact),3));
    
    dEEG.data = fixedEpochsL(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanFixedL = squeeze(mean(dEEG.data(15,:,~isArtifact),3));
    
    allMeanFixed = [meanFixedS' meanFixedM' meanFixedL'];

    dEEG.data = scaledEpochsS(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanScaledS = squeeze(mean(dEEG.data(15,:,~isArtifact),3));
    
    dEEG.data = scaledEpochsM(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanScaledM = squeeze(mean(dEEG.data(15,:,~isArtifact),3));
    
    dEEG.data = scaledEpochsL(:,:,:);
    dEEG.pnts = size(dEEG.data,2); dEEG.trials = size(dEEG.data,3); dEEG = pop_rmbase(dEEG,[],bl);
    isArtifact = find_artifacts(dEEG, maxMin, level, step, lowest);
    meanScaledL = squeeze(mean(dEEG.data(15,:,~isArtifact),3));
    
    allMeanScaled = [meanScaledS' meanScaledM' meanScaledL'];

    allRPre(iParticipant,:,:) = corr(allMeanPre);
    allRFixed(iParticipant,:,:) = corr(allMeanFixed);
    allRScaled(iParticipant,:,:) = corr(allMeanScaled);

    % Compute r-squared
    % Try to predict 'i' using 'j'
    for i = 1:3
        for j = 1:3
            allR2Pre(iParticipant,i,j) = myrsquare(allMeanPre(:,i),allMeanPre(:,j));
            allR2Fixed(iParticipant,i,j) = myrsquare(allMeanFixed(:,i),allMeanFixed(:,j));
            allR2Scaled(iParticipant,i,j) = myrsquare(allMeanScaled(:,i),allMeanScaled(:,j));
        end
    end

%     R = corr(allMean);
%     allR(iParticipant,:,:) = R;

%     figure();
%     subplot(1,3,1); plot(allMeanPre); title('Pre'); legend('S','M','L');
%     subplot(1,3,2); plot(allMeanFixed); title('Fixed'); legend('S','M','L');
%     subplot(1,3,3); plot(allMeanScaled); title('Scaled'); legend('S','M','L');
% %     plot(mean(mEpochs(15,:,:),3)');
% %     plot(mean(lEpochs(15,:,:),3)');
%     pause();
%     disp('*');
end

save(fullfile(resultsFolder,['scalefactor_' whichTask '.mat']));
return;

%%
meanPre = squeeze(mean(allRPre,1))
meanFixed = squeeze(mean(allRFixed,1))
meanScaled = squeeze(mean(allRScaled,1))

subplot(1,3,1); imagesc(meanPre,[0 1]); title('Mixed');
subplot(1,3,2); imagesc(meanFixed,[0 1]); title('Fixed Only');
subplot(1,3,3); imagesc(meanScaled,[0 1]); title('Scaled Only');

collPre = (allRPre(:,1,2) + allRPre(:,1,3) + allRPre(:,2,3)) ./ 3;
collFixed = (allRFixed(:,1,2) + allRFixed(:,1,3) + allRFixed(:,2,3)) ./ 3;
collScaled = (allRScaled(:,1,2) + allRScaled(:,1,3) + allRScaled(:,2,3)) ./ 3;

figure();
notBoxPlot([collPre collFixed collScaled]);
[h,p] = ttest(collFixed,collScaled)
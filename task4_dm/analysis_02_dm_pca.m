% Run PCA on residual EEG
%
% Project: Temporal Scaling
% Other m-files required: EEGLAB

% Author: Cameron Hassall, Department of psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% July, 2022


close all; clear all; clc; rng(41); % for reproducibility

% Analysis settings: 'full', 'fixed', or 'scaled'
whichAnalysis = 'scaled'; % 'full' EEG, 'fixed' only (scaled regressed out), 'scaled' only (fixed regressed out)

% Set data, raw, and results folders, change as needed
rawFolder = '/Users/chassall/Raw/doi_10.5061_dryad.5vb8h__v1';
dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task4_dm/data';
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';
% rawFolder = 'D:\Raw\doi_10.5061_dryad.5vb8h__v1';
% dataFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\task4_dm\data';
% resultsFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\results';

participants = {'02','03','04','05','06','07','08','10','11','12','13','14','15','16','18','19','20','21'};
nParticipants = length(participants);
judgeTimes = [0.512 0.64 0.8 1 1.25; 1.056 1.32 1.65 2.0625 2.5781; 1.6 2 2.5 3.125 3.9063];
subconditions = {'fast','slow'};
srate = 200;
meanArtifacts = [];
aveChannel = 0;
% pcaChannelStrings = {'FC3','FC1','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CP2','CP4','P3','P1','Pz','P2','P4'};
% pcaChannelStrings = {'C1','FC3', 'FC1', 'FCz', 'C3', 'Cz', 'CP3', 'CP1', 'CPz'};
pcaChannelStrings = {'Cz','FC1', 'FCz', 'FC2', 'C1', 'C2', 'CP1', 'CPz', 'CP2'};
pcaChannelStrings = {'Fp1'	'Fp2'	'F3'	'F4'	'C3'	'C4'	'P3'	'P4'	'O1'	'O2'	'F7'	'F8'	'T7'	'T8'	'P7'	'P8'	'Fz'	'Cz'	'Pz'	'Oz'	'AF1'	'AF2'	'FC1'	'FC2'	'CP1'	'CP2'	'PO1'	'PO2'	'FC5'	'FC6'	'CP5'	'CP6'	'F1'	'F2'	'C1'	'C2'	'P1'	'P2'	'AF5'	'AF6'	'FC3'	'FC4'	'CP3'	'CP4'	'PO5'	'PO6'	'F5'	'F6'	'C5'	'C6'	'P5'	'P6'	'AF7'	'AF8'	'FT7'	'FT8'	'PO7'	'PO8'	'FPz'	'FCz'	'CPz'	'NFPz'};

if aveChannel
    whichChannelIndex = 1;
else
    whichChannelString = 'CP3'; % Channel of interest
    whichChannelString = 'Cz'; % Channel of interest
    whichChannelIndex = find(strcmp(pcaChannelStrings,whichChannelString));
end
whichPC = 2; % PC of interest
nQuantiles = 2; % Number of quantile boundaries
epochWindow = [0 0.8];
quantileNames = {'q1','q2','q3'};
withinNames = [{'1'};{'2'};{'3'}];
rmModel = 'q1-q3 ~ 1';
numPCs = 6;
allBinAverages = {};
explainedByParticipant = [];
allPCAScores = [];
allPCACoeff = [];
allPCAScores2 = [];
allDataForPCA = [];

allMeanRTs = [];
allRTs = [];
for iParticipant = 1:nParticipants
    disp(['participant ' participants{iParticipant}]);
    
    % Load the residuals, which will be in EEGLAB format (called rEEG)
    residualsFolder = fullfile(dataFolder, 'derivatives','eegresidual', ['sub-' participants{iParticipant}]);
    residualsFile = ['sub-' participants{iParticipant} '_task-dm_eegresidual_reref'];
    load(fullfile(residualsFolder,residualsFile),'rEEGWithTS','rEEGWithFixed','beta','toExclude');
    
    switch whichAnalysis
        case 'full'
            rEEG = EEG;
        case 'fixed'
            rEEG = rEEGWithFixed;
        case 'scaled'
            rEEG = rEEGWithTS;
    end
    
    pcaChannelIndices = eeg_chaninds(rEEG,pcaChannelStrings);
    
    eventTimes = [rEEG.event.latency];
    eventRTs = eventTimes(2:2:end) - eventTimes(1:2:end-1);
    
    allMeanRTs(iParticipant) = mean(eventRTs);
    allRTs = [allRTs eventRTs];
    
    % Epoch data and find artifacts
    epochedEEG = pop_epoch(rEEG,{1},epochWindow);
    [isArtifact, isArtifactsCT]  = find_artifacts(epochedEEG);
    meanArtifacts(iParticipant) = mean(isArtifact);
    
    % Make quantiles
    quantiles = quantile(eventRTs,nQuantiles);
    theseCategories = nan(1,length(eventRTs));
    for i = 1:length(eventRTs)
        theseCategories(i) = sum(find(eventRTs(i) <= [quantiles max(eventRTs)],1));
    end
    
    % Make bin averages
    binAverages = {};
    binTotals = {};
    for b = 1:nQuantiles+1
        theseTrials = theseCategories == b;
        binAverages{b} = mean(epochedEEG.data(pcaChannelIndices,:,theseTrials & ~isArtifact),3);
        allBinAverages{iParticipant,b} = mean(epochedEEG.data(pcaChannelIndices,:,theseTrials & ~isArtifact),3);
        binCount{b} = sum(theseTrials);
    end
    
    % Do PCA
    pcaScores = [];
    pcaCoeff = {};
    pcaScore = {};
    
    pcaData = [];
    for b = 1:nQuantiles+1
        if aveChannel
            pcaData = [pcaData; mean(binAverages{b},1)];
        else
            pcaData = [pcaData; binAverages{b}];
        end
    end
    
    % For the PCA on all participant data
    allDataForPCA = [allDataForPCA; pcaData];
    
    % Do the PCA for this participant
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(pcaData,'Centered',false);
    explainedByParticipant(iParticipant,:) = EXPLAINED(1:10);
    
%     subplot(1,2,1);
%     plot(COEFF(:,1:2)); legend('1','2');
    
%     % Flip the sign on the first PC as needed
%         if mean(COEFF(:,2)) > 0
%             COEFF(:,2) = -COEFF(:,2);
%             SCORE(:,2) = -SCORE(:,2);
%         end
%     
%         % Flip the sign on the second PC as needed
%         if mean(COEFF(1:end/2,1)) > mean(COEFF(:,1))
%             COEFF(:,1) = -COEFF(:,1);
%             SCORE(:,1) = -SCORE(:,1);
%         end
    
    if mean(COEFF(:,1)) > 0
        COEFF(:,1) = -COEFF(:,1);
        SCORE(:,1) = -SCORE(:,1);
    end
    if mean(COEFF(1:end/2,2)) > mean(COEFF(:,2))
        COEFF(:,2) = -COEFF(:,2);
        SCORE(:,2) = -SCORE(:,2);
    end
    
%         subplot(1,2,2);
%         plot(COEFF(:,1:2)); legend('1','2');
%         pause();
    
    % Store PCA scores for this condition
    if aveChannel
        pcaScores =  reshape(SCORE(:,whichPC),1,[]);
    else
        pcaScores =  reshape(SCORE(:,whichPC),length(pcaChannelIndices),[]);
    end
    pcaCoeff = COEFF;
    pcaScore = SCORE;
    allPCAScores(iParticipant,:,:,:) = pcaScores;
    allPCACoeff(iParticipant,:,:) = pcaCoeff;
end

% Save PCA results for plotting later
quantileData = squeeze(allPCAScores(:,whichChannelIndex,:));
save(fullfile(resultsFolder, 'dm_pca.mat'),'quantileData','nParticipants','explainedByParticipant');
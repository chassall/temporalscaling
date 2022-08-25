% Check preprocessed EEG from the Prediction task
%
% Other m-files required: eeglab

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% 25-Aug-2022

close all; clear all;
ps = [1:3 5 6 7:12 14 16:21];
participants = num2str(ps,'%0.2d,');
participants(end) = [];
participants = strsplit(participants,',');
nParticipants = length(participants);
srate = 200;
testRange = 10000:20000;
task = 'entrainrhythm1b';
for iParticipant = 1:nParticipants
    
    % Load preprocessed data
    preprocessedFolder = ['./data/derivatives/eegprep/sub-' participants{iParticipant}];
    preprocessedFile = ['sub-' participants{iParticipant} '_task-' task '_eegprep'];
    load(fullfile(preprocessedFolder,preprocessedFile), 'EEG');
    
    % Remove ocular components using the results of ICLabel
    eyeLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Eye'));
    eyeClassifications = EEG.etc.ic_classification.ICLabel.classifications(:,eyeLabel);
    [~,eyeI] = sort(eyeClassifications,'descend');
      
    % Pick some data points and plot effect of IC removal
    leftChannel = eeg_chaninds(EEG,'Fp1');
    rightChannel = eeg_chaninds(EEG,'Fp2');
    figure();
    subplot(1,2,1);
    plot(EEG.data(1,testRange)); hold on;
    subplot(1,2,2);
    plot(EEG.data(2,testRange)); hold on;
    for i = 1:6
        thisEyeI = eyeI(1:i);
        testEEG = pop_subcomp(EEG,thisEyeI,0);
        subplot(1,2,1);
        plot(testEEG.data(leftChannel,testRange)); title(participants{iParticipant});
        subplot(1,2,2);
        plot(testEEG.data(rightChannel,testRange));title(participants{iParticipant});
    end
    legend('Original','-IC1','-IC1-2','-IC1-3','-IC1-4','-IC1-5','-IC1-6');
    pause();
    
end

numICs = [1 1 2 2 1 1 2 2 1 1 1 2 1 2 2 1 1 1 ]; % These are manually set based on the above plots


%%
close all; clear all;
ps = [1:3 5 6 7:12 14 16:21];
participants = num2str(ps,'%0.2d,');
participants(end) = [];
participants = strsplit(participants,',');
nParticipants = length(participants);
srate = 200;
testRange = 10000:20000;
task = 'entrainrhythm1b';

maxMin = 150;
level = 150;
step = 40;
lowest = 0.1;
allMeanArtifacts = [];
for iParticipant = 1:nParticipants
    preprocessedFolder = ['../data/derivatives/eegprep/sub-' participants{iParticipant}];
    preprocessedFile = ['sub-' participants{iParticipant} '_task-' task '_eegprep'];
    load(fullfile(preprocessedFolder,preprocessedFile), 'EEG');
    eyeLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Eye'));
    [~,I] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
    whichOnes = find(I == eyeLabel) % Max out of all possible labels
    EEG = pop_subcomp(EEG,whichOnes,0);
    
    eEEG = pop_epoch(EEG,{71,72,81,82,131,132,141,142,231,234,243,244,31,32,41,42},[-0.2,0.6]);
    [isArtifact, isArtifactsCT] = find_artifacts(eEEG, maxMin, level, step, lowest);
    meanArtifacts = mean(isArtifactsCT,2);
    allMeanArtifacts = [allMeanArtifacts meanArtifacts];
end

%%
close all; clear all;
ps = [1:3 5 6 7:12 14 16:21];
ps = 4;
participants = num2str(ps,'%0.2d,');
participants(end) = [];
participants = strsplit(participants,',');
nParticipants = length(participants);
srate = 200;
testRange = 10000:20000;
task = 'entrainrhythm1b';

maxMin = 150;
level = 150;
step = 40;
amplitudeThreshold = 150;
winms = 2000;
lowest = 0.1;
stepms = 1000;
allMeanArtifacts = [];
for iParticipant = 1:nParticipants
    preprocessedFolder = ['../data/derivatives/eegprep/sub-' participants{iParticipant}];
    preprocessedFile = ['sub-' participants{iParticipant} '_task-' task '_eegprep'];
    load(fullfile(preprocessedFolder,preprocessedFile), 'EEG');
    eyeLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Eye'));
    [~,I] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
    whichOnes = find(I == eyeLabel) % Max out of all possible labels
    EEG = pop_subcomp(EEG,whichOnes,0);
    
    % [winrej chanrej] = basicrap(EEG, chanArray, ampth, windowms, stepms, firstdet, fcutoff, forder, thresholdType, numChanThreshold)
    [WinRej, chanrej] = basicrap(EEG, 1:EEG.nbchan, amplitudeThreshold, winms, stepms, 1,[],[],'peak-to-peak',1);
    allMeanArtifacts(iParticipant,:) = mean(chanrej);

end



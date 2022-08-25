% Load and preprocess raw EEG data in BIDS format
%
% Project: Temporal Scaling
% Other m-files required: EEGLAB, find_artifacts.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% March 2020; Last revision: 11-Jul-2022

participants = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};
nParticipants = length(participants);

dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task1-2_productionperception/data';
    
for iParticipant = 1:nParticipants
    
    rng(41); % for reproducibility
    
    %% Load the raw EEG.
    pString = ['sub-' participants{iParticipant}];
    thisRawFolder = fullfile(dataFolder, pString, 'eeg');
    rawFile = [pString '_task-temporalscaling_eeg.set'];
    EEG = pop_loadset(fullfile(thisRawFolder,rawFile));
    
    %% Downsample to 200 Hz.
    EEG = pop_resample(EEG, 200);
    
    %% Apply a bandpass filter (0.1-20 Hz).
    EEG = pop_eegfiltnew(EEG, 0.1, 20);
    
    %% Notch filter
    EEG = pop_eegfiltnew(EEG, 48, 52,[],1);
    
    %% Uncomment for manual inspection - click on "scroll data" then "OK".
    % EEG = pop_select(EEG);
    
    %% Re-reference to the average of the left and right mastoids.
    EEG = pop_reref(EEG, {'LM','RM'});
    
    %% Remove ocular channels from main dataset
    EEG = pop_select(EEG,'nochannel' ,{'LEOG' 'REOG' 'VEOG'});

    %% Make a copy of the full locations file
    fullLocs = EEG.chanlocs;
    
    %% Isolate some data on which to run the ICA
    icaEEG = pop_epoch(EEG,{11,21,31,41,51,61},[0 3]); % 3 s from fixation
    
    %% Remove bad trials from icaEEG (> 1000 uV change).
    badTrialIndex = find_artifacts(icaEEG,500,500,40,0.1);
    icaEEG = pop_select(icaEEG,'notrial',badTrialIndex);
    
    %% Run ICA and get the results.
    % Possible algorithms: 'binica','fastica','runica'.
    icaEEG = pop_runica(icaEEG,'runica');
    icaact = icaEEG.icaact;
    icachansind = icaEEG.icachansind;
    icasphere = icaEEG.icasphere;
    icasplinefile = icaEEG.icasplinefile;
    icaweights = icaEEG.icaweights;
    icawinv = icaEEG.icawinv;
    
    %% Transfer the ICA results from icaEEG to EEG
    EEG.icaact = icaact;
    EEG.icachansind = icachansind;
    EEG.icasphere = icasphere;
    EEG.icasplinefile = icasplinefile;
    EEG.icaweights = icaweights;
    EEG.icawinv = icawinv;
        
    %% Perform IC rejection using the ICLabel EEGLAB extension.
    EEG = iclabel(EEG, 'default');
    
    %% Interpolate missing channels (uncomment if needed).
    EEG = pop_interp(EEG,fullLocs);
    
    %% Save preprocessed EEG
    thisPreprocessedFolder = fullfile(dataFolder, 'derivatives','eegprep',pString);
    if ~exist(thisPreprocessedFolder,'dir'), mkdir(thisPreprocessedFolder); end
    preprocessedFile = [pString '_task-temporalscaling_eegprep'];
    save(fullfile(thisPreprocessedFolder,preprocessedFile),'EEG');
end

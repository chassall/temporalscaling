% Preprocess EEG from the Prediction task
%
% Other m-files required: eeglab

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% 25-Aug-2022

eeglab(); clear all; close all;

rawDir = '/Volumes/DELTA/ccn_lab_raw/temporalprediction'; % Change as needed
ps = [1:12 14:21];

for pi = 1:length(ps)
   

    thisFile = fullfile(rawDir,['EntrainRhythm1B_S1' num2str(ps(pi),'%0.02i') '.bdf']);
    EEG = pop_biosig(thisFile,'ref',[70,71],'refoptions',{'keepref' 'off'}); % Need to re-reference here
    
    if ps(pi) == 4
        otherFile = fullfile(rawDir,['EntrainRhythm1B_S1' num2str(ps(pi),'%0.02i') '_1.bdf']);
        otherEEG = pop_biosig(otherFile,'ref',[70,71],'refoptions',{'keepref' 'off'}); % Need to re-reference here
        EEG = pop_mergeset(otherEEG,EEG);
    elseif ps(pi) == 12
        otherFile = fullfile(rawDir,['EntrainRhythm1B_S1' num2str(ps(pi),'%0.02i') '_2.bdf']);
        otherEEG = pop_biosig(otherFile,'ref',[70,71],'refoptions',{'keepref' 'off'}); % Need to re-reference here
        EEG = pop_mergeset(EEG,otherEEG);
    end
    
    % Fix marker offset
    for i = 1:length(EEG.event)
        if ischar(EEG.event(i).type) EEG.event(i).type = str2num(EEG.event(i).type); end
        EEG.event(i).type = EEG.event(i).type - bin2dec('1111111100000000');
    end
    
    EEG = pop_select(EEG,'nochannel',{'Nose', 'LHEOG', 'RHEOG', 'RVEOGS', 'RVEOGI', 'LVEOGI', 'Ana1', 'Ana2', 'Ana3', 'Ana4', 'Ana5', 'Ana6', 'Ana7', 'Ana8'});
    
    % Assign locs
    EEG.chanlocs = readlocs('temporalprediction64.locs');
    fullLocs = EEG.chanlocs;
    
    %% Remove channel?
    if ps(pi) == 19
        EEG = pop_select(EEG,'nochannel',9);
    end
    
    %% Downsample to 200 Hz
    EEG = pop_resample(EEG, 200);
    
    %% Bandpass filter (0.1-20 Hz for the regular data, 1-20 Hz for the ICA, see "Makoto Pipeline" https://sccn.ucsd.edu/wiki/Makoto%27s_preprocessing_pipeline)
    EEG = pop_eegfiltnew(EEG, 0.1, 20);
    
    %% Notch filter
    EEG = pop_eegfiltnew(EEG, 48, 52,[],1);
    
    %% Isolate some data on which to run the ICA
    icaEEG = pop_epoch(EEG,{71,72,81,82,131,132,141,142,231,234,243,244,31,32,41,42},[-3 0]); % 3-second epoch prior to warning signal
    
    %% Remove bad data from ICA EEG
    badTrialI = any(icaEEG.data < -500,2) | any(icaEEG.data > 500,2); % +/- 500 uV cutoff, aka pretty bad data
    badTrialI = squeeze(any(badTrialI,1));
    icaEEG = pop_select(icaEEG,'notrial',badTrialI);
    
    %% Run ICA
    icaEEG = pop_runica(icaEEG,'runica'); % Possibilities: 'binica','fastica','runica'
    
    %% Transfer the ICA results from icaEEG to EEG
    preICA = EEG; % For later comparison
    EEG.icaact = icaEEG.icaact;
    EEG.icachansind = icaEEG.icachansind;
    EEG.icasphere = icaEEG.icasphere;
    EEG.icasplinefile = icaEEG.icasplinefile;
    EEG.icaweights = icaEEG.icaweights;
    EEG.icawinv = icaEEG.icawinv;
    
    %% Perform IC rejection using ICLabel scores and r.v. from dipole fitting.
    
    EEG = iclabel(EEG, 'default');
    eyeLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Eye'));
    eyeI  = find(EEG.etc.ic_classification.ICLabel.classifications(:,eyeLabel) >= 0.7);
    
    %% Remove ocular component
    % EEG = pop_subcomp(EEG,eyeI,0);
    
    %% Interpolate missing channels
    EEG = pop_interp(EEG,fullLocs);
   
    %% Save
    task = 'entrainrhythm1b';
    preprocessedFolder = ['./data/derivatives/eegprep/sub-' num2str(ps(pi),'%0.02i')];
    preprocessedFile = ['sub-' num2str(ps(pi),'%0.02i') '_task-' task '_eegprep.mat'];
    preFile = fullfile(preprocessedFolder,preprocessedFile);
    save(preFile,'EEG');
end
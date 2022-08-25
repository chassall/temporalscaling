% Load source EEG data in CURRY format and convert to BIDS
%
% Project: Temporal Scaling
% Other m-files required: EEGLAB with loadcurry and bids-matlab-tools extensions
% bids-matlab toolbox (https://github.com/bids-standard/bids-matlab)
% Other files required: neuroscan37.locs

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% June 2020; Last revision: 11-Jul-2022

% Folders
% 1. data_raw: raw EEG and behavioural, not in BIDS format (although the
% EEG is in a BIDS-friendly folder called sourcedata)
% 2. data_temp: a temporary folder where loaded raw EEG will be saved as
% BIDS-friendly .set files
% 3. data: the BIDS folder containing EEG, behavioural, and (eventually)
% derivatives (this will be uploaded to OpenNeuro)

% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};

% Folder locations (change as needed)
rawFolder = '/Users/chassall/Library/CloudStorage/OneDrive-Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task1-2_productionperception/data_raw';
tempFolder = '/Users/chassall/Library/CloudStorage/OneDrive-Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task1-2_productionperception/data_temp';
bidsFolder = '/Users/chassall/Library/CloudStorage/OneDrive-Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task1-2_productionperception/data';

% Make the temp folder if it doesn't exist
if ~exist(tempFolder,'dir')
    mkdir(tempFolder);
end

% Loop through participants, load their raw data, save as .set in the temp
% folder
for iParticipant = 1:length(ps)
    
    % Load the CURRY source data
    thisSourceFolder = [rawFolder '/sourcedata/sub-' ps{iParticipant} '/eeg' ];
    sourceFile = ['sub-' ps{iParticipant} '_task-temporalscaling_eeg.cdt'];
    EEG = loadcurry(fullfile(thisSourceFolder, sourceFile));
    
    % Remove trigger channel
    EEG.data(end,:) = [];
    EEG.nbchan = EEG.nbchan - 1;
    
    % Add in the reference channel.
    EEG.data = [EEG.data; zeros(1,EEG.pnts)];
    EEG.nbchan = EEG.nbchan + 1;
    
    % Load/set channel locations file.
    EEG.chanlocs = readlocs('neuroscan37.locs');
    
    % Save the data in EEGLAB format.
    thisRawFolder = [tempFolder '/sub-' ps{iParticipant} '/eeg' ];
    if ~exist(thisRawFolder,'dir'), mkdir(thisRawFolder); end
    rawFile = ['sub-' ps{iParticipant} '_task-temporalscaling_eeg.set'];
    pop_saveset(EEG,'filename',fullfile(thisRawFolder, rawFile));
end

%% Make and save all necessary BIDS files using bids_export()

% Make the BIDS folder if it doesn't exist
if ~exist(bidsFolder,'dir')
    mkdir(bidsFolder);
end

% Set data file names 
for p = 1:length(participants)
    data(p).file = {fullfile(tempFolder, ['sub-' ps{p}], 'eeg', ['sub-' ps{p} '_task-temporalscaling_eeg.set'])};
end

% dataset_description.json
generalInfo.Name = 'Temporal Scaling';
generalInfo.Authors = {'Cameron D. Hassall', 'Laurence T. Hunt'};
generalInfo.License = 'CC0';
generalInfo.DatasetType = 'raw';
generalInfo.DatasetDOI = '';

% participants.tsv 
pInfo = {'participant_id', 'datetime', 'age', 'sex', 'handedness';
    'sub-01', '20-Feb-2020 14:55:06', '25', 'M',  'R';
    'sub-02', '27-Jan-2020 12:13:09', '21', 'F',  'R';
    'sub-03', '01-Feb-2020 16:33:09', '24' , 'M',  'R';
    'sub-04', '04-Feb-2020 10:53:23', '21', 'F',  'R';
    'sub-05', '05-Feb-2020 10:15:49', '25', 'F',  'L';
    'sub-06', '05-Feb-2020 11:57:35', '20', 'M',  'R';
    'sub-07', '05-Feb-2020 14:57:18', '21', 'M',  'R';
    'sub-08', '12-Feb-2020 15:25:55', '30', 'M',  'L';
    'sub-09', '15-Feb-2020 14:58:24', '23', 'F',  'R';
    'sub-10', '18-Feb-2020 10:13:32', '24', 'F',  'R';
    'sub-11', '05-May-2022 11:28:12', '64', 'F',  'L';
    'sub-12', '06-May-2022 11:21:12', '65', 'F',  'R';
    'sub-13', '11-May-2022 11:25:21', '40', 'F',  'L';
    'sub-14', '11-May-2022 14:10:05', '67', 'F',  'R';
    'sub-15', '12-May-2022 12:05:38', '71', 'M',  'R';
    'sub-16', '12-May-2022 16:03:13', '77', 'M',  'R';
    'sub-17', '18-May-2022 12:07:47', '66', 'F',  'L';
    'sub-18', '26-May-2022 15:23:19', '53', 'F',  'LR';
    'sub-19', '01-Jun-2022 11:05:26', '49', 'F',  'R';
    'sub-20', '24-Jun-2022 10:33:12', '39', 'F',  'R'};

% participants.json
pInfoDesc.participant_id.Description = 'participant number';
pInfoDesc.datetime.Description = 'date and time at start of task';
pInfoDesc.age.Description = 'self-reported age of participant';
pInfoDesc.sex.Description = 'self-reported sex of participant';
pInfoDesc.sex.Levels.M = 'male';
pInfoDesc.sex.Levels.F = 'female';
pInfoDesc.handedness.Description = 'self-reported handedness of participant';
pInfoDesc.handedness.Levels.L = 'left-handed';
pInfoDesc.handedness.Levels.R = 'right-handed';
pInfoDesc.handedness.Levels.LR = 'ambidextrous';

% task-temporalscaling_events.json
eInfoDesc.type.Description = 'Event value (int)';
eInfoDesc.type.Levels.x2 = 'Pre-block metronome beep';
eInfoDesc.type.Levels.x11 = 'Fixation cross (production task, fast tempo)';
eInfoDesc.type.Levels.x12 = 'Start of trial beep (production task, fast tempo)';
eInfoDesc.type.Levels.x13 = 'Participant response (production task, fast tempo)';
eInfoDesc.type.Levels.x14 = 'Correct feedback (production task, fast tempo)';
eInfoDesc.type.Levels.x15 = 'Early feedback (production task, fast tempo)';
eInfoDesc.type.Levels.x16 = 'Late feedback (production task, fast tempo)';
eInfoDesc.type.Levels.x21 = 'Fixation cross (production task, medium tempo)';
eInfoDesc.type.Levels.x22 = 'Start of trial beep (production task, medium tempo)';
eInfoDesc.type.Levels.x23 = 'Participant response (production task, medium tempo)';
eInfoDesc.type.Levels.x24 = 'Correct feedback (production task, medium tempo)';
eInfoDesc.type.Levels.x25 = 'Early feedback (production task, medium tempo)';
eInfoDesc.type.Levels.x26 = 'Late feedback (production task, medium tempo)';
eInfoDesc.type.Levels.x31 = 'Fixation cross (production task, slow tempo)';
eInfoDesc.type.Levels.x32 = 'Start of trial beep (production task, slow tempo)';
eInfoDesc.type.Levels.x33 = 'Participant response (production task, slow tempo)';
eInfoDesc.type.Levels.x34 = 'Correct feedback (production task, slow tempo)';
eInfoDesc.type.Levels.x35 = 'Early feedback (production task, slow tempo)';
eInfoDesc.type.Levels.x36 = 'Late feedback (production task, slow tempo)';
eInfoDesc.type.Levels.x41 = 'Fixation cross (perception task, fast tempo)';
eInfoDesc.type.Levels.x42 = 'Start of trial beep (perception task, fast tempo)';
eInfoDesc.type.Levels.x43 = 'Participant response (perception task, fast tempo)';
eInfoDesc.type.Levels.x44 = 'Correct feedback (perception task, fast tempo)';
eInfoDesc.type.Levels.x45 = 'Early feedback (perception task, fast tempo)';
eInfoDesc.type.Levels.x46 = 'Late feedback (perception task, fast tempo)';
eInfoDesc.type.Levels.x51 = 'Fixation cross (perception task, medium tempo)';
eInfoDesc.type.Levels.x52 = 'Start of trial beep (perception task, medium tempo)';
eInfoDesc.type.Levels.x53 = 'Participant response (perception task, medium tempo)';
eInfoDesc.type.Levels.x54 = 'Correct feedback (perception task, medium tempo)';
eInfoDesc.type.Levels.x55 = 'Early feedback (perception task, medium tempo)';
eInfoDesc.type.Levels.x56 = 'Late feedback (perception task, medium tempo)';
eInfoDesc.type.Levels.x61 = 'Fixation cross (perception task, slow tempo)';
eInfoDesc.type.Levels.x62 = 'Start of trial beep (perception task, slow tempo)';
eInfoDesc.type.Levels.x63 = 'Participant response (perception task, slow tempo)';
eInfoDesc.type.Levels.x64 = 'Correct feedback (perception task, slow tempo)';
eInfoDesc.type.Levels.x65 = 'Early feedback (perception task, slow tempo)';
eInfoDesc.type.Levels.x66 = 'Late feedback (perception task, slow tempo)';
eInfoDesc.latency.Description = 'Event onset';
eInfoDesc.latency.Units = 'samples';
eInfoDesc.urevent.Description = 'Event number';

% README
README = sprintf('# Temporal Scaling\n\nTwenty participants learned three temporal intervals. There were two subtasks, randomly interleaved. In the Production task, participants produced either a short, medium, or long temporal interval. In the Perception task, participants judged a computer-produced interval as correct or incorrect (again, for a short, medium, or long temporal interval). In both tasks participants received visual feedback (a checkmark or x).\n\nPreprint: https://doi.org/10.1101/2020.12.11.421180');

% We won't do a CHANGES file because this will be generated by OpenNeuro

% sub-XX_task-temporalscaling_eeg.json
tInfo.InstitutionAddress = 'Warneford Hospital, Warneford Lane';
tInfo.InstitutionName = 'University of Oxford';
tInfo.InstitutionalDepartmentName = 'Department of Psychiatry';
tInfo.PowerLineFrequency = 50;
tInfo.ManufacturersModelName = 'actiCHamp Plus';

% sub-XX_task-temporalscaling_channels.tsv
chanlocs = '/Users/chassall/Library/CloudStorage/OneDrive-Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task1-2_productionperception/neuroscan37.locs';

% Run bids_export
bids_export(data, 'targetdir', bidsFolder, 'taskName', 'temporalscaling', 'gInfo', generalInfo, 'pInfo', pInfo, 'pInfoDesc', pInfoDesc, 'eInfoDesc', eInfoDesc, 'README', README,'tInfo', tInfo, 'chanlocs', chanlocs);

%% Load the behavioural files and save as TSVs
for p = 1:length(ps)

    % Load behavioural data from rawFolder/beh
    pBIDSString = ['sub-' ps{p}];
    thisFile = fullfile(rawFolder, 'beh',[pBIDSString '_task-temporalscaling_beh.txt']);
    thisData = load(thisFile);

    % Make a struct out of the behavioural data
    beh.block = thisData(:,1);
    beh.trial = thisData(:,2);
    beh.preBeepTime = thisData(:,3); 
    beh.responseTime = thisData(:,4); 
    beh.computerRT = thisData(:,5); 
    beh.computerRespCondition = thisData(:,6); 
    beh.judgementMapping = thisData(:,7); 
    beh.participantResponse = thisData(:,8); 
    beh.participantJudgement = thisData(:,9); 
    beh.guessTime = thisData(:,10); 
    beh.preFeedbackTime = thisData(:,11); 
    beh.fbCondition = thisData(:,12); 
    beh.totalPoints = thisData(:,13); 
    beh.thisTrialMargin = thisData(:,14);

    % Make /beh folder for this participant
    behFolder = fullfile(bidsFolder,pBIDSString,'beh');
    if ~exist(behFolder)
        mkdir(behFolder);
    end

    behFile = [pBIDSString '_task-temporalscaling_beh.tsv'];
    bids.util.tsvwrite(fullfile(behFolder,behFile), beh);
end

%% write beh json for each participant
for p = 1:length(ps)

    pString = ['sub-' ps{p}];
    behJSONFile = fullfile(bidsFolder, pString,'beh',[pString '_task-temporalscaling_beh.json']);
    
    bInfoDesc.block.Description = 'Block type (integer)';
    bInfoDesc.block.Levels.x1 = 'Production task, fast tempo';
    bInfoDesc.block.Levels.x2 = 'Production task, medium tempo';
    bInfoDesc.block.Levels.x3 = 'Production task, slow tempo';
    bInfoDesc.block.Levels.x4 = 'Perception task, fast tempo';
    bInfoDesc.block.Levels.x5 = 'Perception task, medium tempo';
    bInfoDesc.block.Levels.x6 = 'Perception task, slow tempo';
    bInfoDesc.trial.Description = 'Trial number (integer)';
    bInfoDesc.preBeepTime.Description = 'Fixation time before cue onset (float)';
    bInfoDesc.preBeepTime.Units = 'seconds';
    bInfoDesc.responseTime.Description = 'Production task: Response time (float)';
    bInfoDesc.responseTime.Units = 'seconds';
    bInfoDesc.computerRT.Description = 'Perception task: Computer-produced interval (float)';
    bInfoDesc.computerRT.Units = 'seconds';
    bInfoDesc.computerRespCondition.Description = 'Perception task: Condition number (int)';
    bInfoDesc.computerRespCondition.Levels.x1 = 'Very early';
    bInfoDesc.computerRespCondition.Levels.x2 = 'Early';
    bInfoDesc.computerRespCondition.Levels.x3 = 'On time';
    bInfoDesc.computerRespCondition.Levels.x4 = 'Late';
    bInfoDesc.computerRespCondition.Levels.x5 = 'Very late';
    bInfoDesc.judgementMapping.Description = 'Perception task: Mapping from keypress to yes/no judgement (int)';
    bInfoDesc.judgementMapping.Levels.x1 = 'left = yes, right = no';
    bInfoDesc.judgementMapping.Levels.x2 = 'left = no, right = yes';
    bInfoDesc.participantResponse.Description = 'Perception task: Participant left/right response (int)';
    bInfoDesc.participantResponse.Levels.x1 = 'left';
    bInfoDesc.participantResponse.Levels.x2 = 'right';
    bInfoDesc.participantJudgement.Description = 'Perception task: Participant yes/no judgement (int)';
    bInfoDesc.participantJudgement.Levels.x1 = 'yes';
    bInfoDesc.participantJudgement.Levels.x2 = 'no';
    bInfoDesc.guessTime.Description = 'Perception task: Participant judgment time';
    bInfoDesc.guessTime.Units = 'seconds';
    bInfoDesc.preFeedbackTime.Description = 'Fixation time before feedback onset (float)';
    bInfoDesc.preFeedbackTime.Units = 'seconds';
    bInfoDesc.fbCondition.Description = 'Feedback type (int)';
    bInfoDesc.fbCondition.Levels.x1 = 'Perception task: Incorrect judgement';
    bInfoDesc.fbCondition.Levels.x2 = 'Perception task: Correct judgement';
    bInfoDesc.fbCondition.Levels.x3 = 'Production task: Early response';
    bInfoDesc.fbCondition.Levels.x4 = 'Production task: On time response';
    bInfoDesc.fbCondition.Levels.x5 = 'Production task: Late response';
    bInfoDesc.totalPoints.Description = 'Point total (int)';
    bInfoDesc.thisTrialMargin.Description = 'Production task: timing margin (float)';
    bInfoDesc.thisTrialMargin.Units = 'seconds';
    options.indent = '  '; % Adds white space, easier to read
    bids.util.jsonwrite(behJSONFile,bInfoDesc,options);
end

%% Manual steps after running this script
% - delete data_temp
% - delete /data/stimuli
% - create /data/code folder and copy into it this script and task script

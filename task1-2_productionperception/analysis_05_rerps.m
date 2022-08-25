% Load preprocessed EEG and make regression ERPs (rERPs)
%
% Project: Temporal Scaling
% Other m-files required: EEGLAB, basicrap.m, pinv_reg.m

% Author: Cameron Hassall, Department of psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% March 2020; Last revision: 12-Jul-2022

close all; clear all; clc; rng(41); % for reproducibility

% Set data and results folders (change as needed)
dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task1-2_productionperception/data';
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

% Task flag ('production' or 'perception')
whichTask = 'perception';
probeTimes = [0.512 0.64 0.8 1 1.25; 1.056 1.32 1.65 2.0625 2.5781; 1.6 2 2.5 3.125 3.9063];
subconditions = {'very early','early','on time','late','very late'};
subConditionsToInclude = {'very early','early','on time','late','very late'};

% Participant information
participants = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};
nParticipants = length(participants);

% EEG parameters
doLaplacian = 0;
srate = 200;
amplitudeThreshold = 150; % Artifact amplitude threshold
winms = 2000; % Artifact window size
stepms = 1000; % Artifact detection step size

% Behavioural parameters
behavLimits = [0.2 5]; % For defining behavioural artifacts

% Method for stretching/compressing the scaled regressors
stretchMethod = 'box'; % 'box', 'triangle', 'cubic', or 'nearest'

% Regularization parameters
doReg = 1;
loadPrevReg = 0; % Re-use previous regularization result (saves time)
k = 10; % Number of fold for cross-validation
regChannels = {'FCz'}; % Do the cross-validation at these channels only
regtype = 'onediff'; % Use a finite difference approximation of the first derivative
lambdas = [0.001 0.01 0 1 10 100 1000 10000 100000 1000000 10000000]; % Reg. parameters to try
allCVErrors = nan(length(participants),32,length(lambdas)); % To store cross-validation errors

allVIF = []; % Variance inflation factor
allX = {}; % Design matrices

% Task-specific parameters
switch whichTask
    case 'production'
        
        % Task details
        conditionNumbers = [1 2 3]; % Block codes in behavioural file
        
        % Stim-locked component
        stimTimelimits = [-0.2,0.8]; % Time window
        stimScalePntLimits = srate * stimTimelimits;
        stimScalePnts = srate * stimTimelimits(1):srate *stimTimelimits(2);
        numStimPnts = length(stimScalePnts);
        stimBL = [-0.2 0];
        
        % Response-locked component
        respTimelimits = [-0.8,0.2]; % Time window
        % respTimelimits = [-3,0.6]; % Time window
        respScalePntLimits = srate * respTimelimits;
        respScalePnts = srate * respTimelimits(1): srate *respTimelimits(2);
        numRespPnts = length(respScalePnts);
        respBL = [-0.8 -0.750]; % Baseline, in seconds
        
        % Time-scaled component
        scaleLength = 330; % In data points
            
        % GLM model type
        % 1 = scaled component spans beep to beep
        % 2 = scaled component spans beep to target time
        % 3 = scaled component spans beep to target time, but gets "cut
        % off" if the beep is early (in line with CNV literature)
        whichGLMModel = 1;
        
    case 'perception'
        
        % Task details
        conditionNumbers = [4 5 6]; % Block codes in behavioural file
        
        % Stim-locked component
        stimTimelimits = [-0.2,0.8]; % Stim start/end times
        stimScalePntLimits = srate * stimTimelimits;
        stimScalePnts = srate * stimTimelimits(1):srate *stimTimelimits(2); % Stim start/end point
        numStimPnts = length(stimScalePnts);
        stimBL = [-0.2 0]; % Baseline, in seconds
        
        % Probe-locked component
        respTimelimits = [-0.8,0.2]; % Response start/end times
        respScalePntLimits = srate * respTimelimits;
        respScalePnts = srate * respTimelimits(1): srate *respTimelimits(2); % Response start/end points
        numRespPnts = length(respScalePnts);
        respBL = [-0.6 -0.550]; % Baseline, in seconds
        
        % Time-scaled component
        scaleLength = 330; % Length of temporally-scaled (TS) component
        
        % GLM model type
        % 1 = scaled component spans beep to beep
        % 2 = scaled component spans beep to target time
        % 3 = scaled component spans beep to target time, but gets "cut
        % off" (in line with CNV literature)
        whichGLMModel = 3;
end
        
% Load previous regularization results? (saves having to do it again)
if loadPrevReg
    load(fullfile(resultsFolder,['rerp_results_' whichTask '_task_' stretchMethod '_' num2str(doLaplacian) '.mat']));
    % load(['../results/'  'rerp_results_' whichTask '_task_' stretchMethod '_' num2str(doLaplacian) '.mat']);
end

tsBreakpoints = [numStimPnts (numStimPnts+numRespPnts)]; % Since the betas will be appended, we need to know where the boundaries are

%% Output variables

% Betas
allBetasStim = nan(length(participants),32,numStimPnts);
allBetasResp = nan(length(participants),32,numRespPnts);
allBetasTS  = nan(length(participants),32,scaleLength);
allRSquared = nan(length(participants),32);
allRSquaredFixed = nan(length(participants),32);

% Residuals
allResiduals = [];
allStimRs = [];
allRespRs = [];

% Artifacts
allContinuousEEGArtifacts = nan(length(participants),1);
allContinuousBehaviouralArtifacts = nan(length(participants),1);

% Mean proportion of artifacts of either type
allContinuousArtifacts = nan(length(participants),1);

%% Participant Loop
% Loop through all participants, run the GLM, save the results
for iParticipant = 1:nParticipants
    
    % Load preprocessed data
    subString = ['sub-' participants{iParticipant}];
    % preprocessedFolder = ['./data/derivatives/eegprep/sub-' participants{iParticipant}];
    preprocessedFolder = fullfile(dataFolder,'derivatives','eegprep',subString);
    preprocessedFile = [subString '_task-temporalscaling_eegprep'];
    disp(['loading ' preprocessedFile]);
    load(fullfile(preprocessedFolder,preprocessedFile), 'EEG');
    
    % Round latencies, as some may be non-integers due to resampling
    for i = 1:length(EEG.event)
        EEG.event(i).latency = round(EEG.event(i).latency);
    end
    
    iRegChannels = eeg_chaninds(EEG,regChannels);
    
    % Remove ocular components using the results of ICLabel
    eyeLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Eye'));
    [~,I] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
    whichOnes = find(I == eyeLabel); % Max out of all possible labels
    EEG = pop_subcomp(EEG,whichOnes,0);
    
    % Make a copy of EEG that will eventually contain the continuous residual
    rEEGWithTS = EEG; % scaled component only (no fixed)
    rEEGWithFixed = EEG;
    
    %  Load behavioural data and get RTs (only needed in the "judge" task)
    % behFolder = ['./data/rawdata/sub-' participants{iParticipant} '/beh'];
    behFolder = fullfile(dataFolder,subString,'beh');
    behFile = [subString '_task-temporalscaling_beh.tsv'];
    opts = delimitedTextImportOptions("NumVariables", 14);
    opts.Delimiter = "\t";
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    thisData = readtable(fullfile(behFolder,behFile), opts);
    thisData = table2array(thisData); % Do this because we originally worked with an array

    % thisData = load(fullfile(behFolder,behFile));
    judgeData = thisData(thisData(:,1)>=4,:);
    
    % Add string labels to EEG data structure
    perceptionTrialCount = 1;
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
                thisRT = judgeData(perceptionTrialCount,5);
                EEG.event(iEvent).subcondition = subconditions{find(thisRT==probeTimes(1,:))};
                perceptionTrialCount = perceptionTrialCount + 1;
                EEG.event(iEvent).task = 'perception';
                EEG.event(iEvent).type = 'response';
                EEG.event(iEvent).condition = 'short';
            case 52
                EEG.event(iEvent).task = 'perception';
                EEG.event(iEvent).type = 'beep';
                EEG.event(iEvent).condition = 'medium';
            case 53
                thisRT = judgeData(perceptionTrialCount,5);
                EEG.event(iEvent).subcondition = subconditions{find(thisRT==probeTimes(2,:))};
                perceptionTrialCount = perceptionTrialCount + 1;
                EEG.event(iEvent).task = 'perception';
                EEG.event(iEvent).type = 'response';
                EEG.event(iEvent).condition = 'medium';
            case 62
                EEG.event(iEvent).task = 'perception';
                EEG.event(iEvent).type = 'beep';
                EEG.event(iEvent).condition = 'long';
            case 63
                thisRT = judgeData(perceptionTrialCount,5);
                EEG.event(iEvent).subcondition = subconditions{find(thisRT==probeTimes(3,:))};
                perceptionTrialCount = perceptionTrialCount + 1;
                EEG.event(iEvent).task = 'perception';
                EEG.event(iEvent).type = 'response';
                EEG.event(iEvent).condition = 'long';
            otherwise % No change
                EEG.event(iEvent).type = num2str(EEG.event(iEvent).type);
        end
    end
    
    % Artifact detection on the continuous EEG
    % This will identify windows of width winms that contain artifacts
    [WinRej, chanrej] = basicrap(EEG, 1:EEG.nbchan, amplitudeThreshold, winms, stepms);
    
    % Flag EEG artifacts
    toExclude = [];
    isEEGArtifact = false(1,size(EEG.data,2));
    for a = 1:size(WinRej,1)
        toExclude = [toExclude WinRej(a,1):WinRej(a,2)];
    end
    isEEGArtifact(toExclude) = true;
    
    % Also flag samples that are behavioural artifacts
    bArtifactTrials = [];
    bArtifactPoints = [];
    trialCount = 1;
    for iEvent = 1:length(EEG.event)
        if strcmp(EEG.event(iEvent).task,whichTask) && strcmp(EEG.event(iEvent).type,'beep')
            thisBeepLatency = EEG.event(iEvent).latency;
            thisResponseLatency = EEG.event(iEvent+1).latency;
            thisRT = (thisResponseLatency-thisBeepLatency)*(1/srate);
            isArtifact = thisRT < behavLimits(1) || thisRT > behavLimits(2);
            if isArtifact
                artifactStart = thisBeepLatency + stimScalePntLimits(1);
                artifactEnd = thisResponseLatency + respScalePntLimits(2);
                bArtifactPoints = [bArtifactPoints artifactStart:artifactEnd];
                bArtifactTrials = [bArtifactTrials trialCount];
            end
            trialCount = trialCount + 1;
        end
    end
    isBehaviouralArtifact = false(1,size(EEG.data,2));
    isBehaviouralArtifact(bArtifactPoints) = true;
    
    % Make the design matrix X: [STIM RESPONSE SCALED]
    
    % Fixed-time components (cue and response/probe)
    stimX = sparse(EEG.pnts,numStimPnts);
    respX = sparse(EEG.pnts,numRespPnts);
    trialCount = 1;
    for iEvent = 1:length(EEG.event)
        thisLatency = EEG.event(iEvent).latency;
        % Is this an event of interest?
        if strcmp(EEG.event(iEvent).task,whichTask) && strcmp(EEG.event(iEvent).type,'beep')
            for j = 1:numStimPnts
                stimX(thisLatency+stimScalePnts(j),j) = 1;
            end
        elseif strcmp(EEG.event(iEvent).task,whichTask) && strcmp(EEG.event(iEvent).type,'response')
            for j = 1:numRespPnts
                respX(thisLatency+respScalePnts(j),j) = 1;
            end
            trialCount = trialCount + 1;
        end
    end
    
    % Scaled-time component
    scaleX = sparse(EEG.pnts,scaleLength);
    scaleEvents = {'beep','response'};
    scaleStart = NaN;
    scaleEnd = NaN;
    trialCount = 1;
    for iEvent = 1:length(EEG.event)
        if strcmp(EEG.event(iEvent).type,scaleEvents{1}) && strcmp(EEG.event(iEvent).task,whichTask)
            scaleStart = EEG.event(iEvent).latency;
        elseif strcmp(EEG.event(iEvent).type,scaleEvents{2}) && strcmp(EEG.event(iEvent).task,whichTask)
            % Find the target time, in data points
            thisTarget = NaN;
            if strcmp(EEG.event(iEvent).condition,'short')
                thisTarget = 0.8*srate;
            elseif strcmp(EEG.event(iEvent).condition,'medium')
                thisTarget = 1.65*srate;
            elseif strcmp(EEG.event(iEvent).condition,'long')
                thisTarget = 2.5*srate;
            else
                error('unrecognized condition');
            end
            % Set the endpoint for the scaling, depending on the model
            switch whichGLMModel
                case 1
                    scaleEnd = EEG.event(iEvent).latency;
                case {2,3}
                    scaleEnd = scaleStart + thisTarget;
                    computerResponseTime = EEG.event(iEvent).latency - scaleStart;
            end
            
            % Start/end points
            whichPnts = scaleStart:scaleEnd;
            
            % Scaled-time regressors. Stretch/compress an identity matrix.
            thisBasis = imresize(eye(scaleLength),[length(whichPnts),scaleLength],stretchMethod);
            
            % Model 3: truncate TS component at the time the second beep occurs
            if whichGLMModel == 3
                thisBasis(computerResponseTime+1:end,:) = 0;
            end
            scaleX(whichPnts,:) = thisBasis;
            trialCount = trialCount + 1; % Increment trial counter
            scaleStart = NaN; % Done with this trial - set this to NaN to avoid mixups
            
        end
    end
    
    % Construct design matrix
    XNoTS = [stimX respX]; % Fixed-time components only (no time-scaling)
    XNoFixed = [scaleX];
    X = [stimX respX scaleX]; % Fixed and scaled components
    XFixed = [stimX respX];
    hasFixed = any(XFixed,2);
    if doLaplacian
        [data,G,H] = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    else
        data = EEG.data;
    end
    
    % Store an copy for the purpose of constucting the residuals
    oX = X; % Original X
    oData = data; % Original EEG
    
    allX{iParticipant} = X; % Store design matrix for each participant (new for PNAS submission)
    
    % Flag empty rows
    isEmptyRow = full(sum(abs(X),2) == 0)';
    
    % Check to see what proportion of included samples are artifacts
    isTrialEEGArtifact = isEEGArtifact & ~isEmptyRow;
    allContinuousEEGArtifacts(iParticipant) = mean(double(isTrialEEGArtifact));
    isTrialBehavArtifact = isBehaviouralArtifact & ~isEmptyRow;
    allContinuousBehaviouralArtifacts(iParticipant) = mean(double(isTrialBehavArtifact));
    
    isEitherArtifact = isTrialEEGArtifact | isTrialBehavArtifact;
    isEitherArtifact(isEmptyRow) = [];
    allContinuousArtifacts(iParticipant) = mean(isEitherArtifact);
    
    % Exclude rows that are empty, an EEG artifact, or a behavioural
    % artifact
    toExclude = isEmptyRow | isTrialEEGArtifact | isTrialBehavArtifact;
    preX = X;
    preX(isEmptyRow,:) = [];
    X(toExclude,:) = [];
    XFixed(toExclude,:) = [];
    hasFixed(toExclude) = [];
    data(:,toExclude) = [];
    fullData = data;
    dataNoTS = data;
    fullX = X;
    beta = nan(size(X,2),EEG.nbchan);
    betaFixed = nan(size(XNoTS,2),EEG.nbchan);

    if doReg
        
        if ~loadPrevReg
            % Tikhonov regularization
            allErrors = nan(EEG.nbchan,length(lambdas));
            tempBetas = [];
            for li = 1:length(lambdas)
                lambda = lambdas(li);
                cv = cvpartition(size(data,2),'KFold',k);
                theseErrors = nan(1,k);
                % Channel loop
                for c = iRegChannels
                    % Cross-validation loop
                    for ki = 1:k
                        % Train
                        thisTrainingI = training(cv,ki);
                        thisTestingI = test(cv,ki);
                        thisPinv = pinv_reg(X(thisTrainingI,:),lambda,regtype,tsBreakpoints);
                        thisBeta = thisPinv * data(c,thisTrainingI)';
                        % Test
                        thisResidual = data(c,thisTestingI)' - full(X(thisTestingI,:)) * thisBeta;
                        thisError = sum(thisResidual.*thisResidual);
                        theseErrors(ki) = thisError;
                    end
                    allErrors(c,li) = mean(theseErrors);
                end
            end
            allCVErrors(iParticipant,:,:) = allErrors;
        end
        
        % Examine errors to pick the best lambda
        [M,I] = min(allCVErrors(iParticipant,iRegChannels,:),[],3); % Might break for more reg channels
        bestLambda = lambdas(I);
        
        % Do regression
        pDM = pinv_reg(X,bestLambda,regtype,tsBreakpoints);
        pDMFixed = pinv_reg(XFixed,bestLambda,regtype,tsBreakpoints(1));
        for c = 1:EEG.nbchan
            beta(:,c) = pDM * data(c,:)';
            betaFixed(:,c) = pDMFixed * data(c,:)';
        end
        
    else
        
        lsmriterations = 400;
        for c = 1:EEG.nbchan
            [beta(:,c),~,~] = lsmr(X,double(data(c,:)'),[],10^-8,10^-8,[],lsmriterations);
        end
    end
    
    % Compute error
    residual = data - (X * beta)';
    residual(:,~hasFixed) = []; % Only count samples that are represented in both analyses (to make it fair)
    rsquared = sum(residual.* residual,2);
    residualFixed = data - (XFixed*betaFixed)';
    residualFixed(:,~hasFixed) = []; 
    rquaredFixed = sum(residualFixed.* residualFixed,2);
    allRSquared(iParticipant,:) = rsquared;
    allRSquaredFixed(iParticipant,:) = rquaredFixed;
    
    % Compute residuals (to be stored)
    residualsWithTS = oData - (XNoTS * beta(1:size(XNoTS,2),:))';
    rEEGWithTS.data = residualsWithTS; % TS only
    
    residualsWithFixed = oData - (XNoFixed * beta(size(XNoTS,2)+1:end,:))';
    rEEGWithFixed.data = residualsWithFixed; % Fixed only
    
    % Store betas
    stimBeta = beta(1:numStimPnts,:)';
    respBeta = beta(numStimPnts+1:numStimPnts+numRespPnts,:)';
    tsBeta = beta(numStimPnts+numRespPnts+1:numStimPnts+numRespPnts+scaleLength,:)';
    allBetasStim(iParticipant,:,:) = stimBeta;
    allBetasResp(iParticipant,:,:) = respBeta;
    allBetasTS(iParticipant,:,:) = tsBeta;
    
    % Save the residuals for this participant
    if doReg || loadPrevReg
        if doLaplacian
            % residualsFolder = ['./data/derivatives/eegresiduallap/sub-' participants{iParticipant}];
            residualsFolder = fullfile(dataFolder,'derivatives/eegresiduallap',subString);
        else
            % residualsFolder = ['./data/derivatives/eegresidual/sub-' participants{iParticipant}];
            residualsFolder = fullfile(dataFolder,'derivatives/eegresidual',subString);
        end
        if ~exist(residualsFolder,'dir'), mkdir(residualsFolder); end
        residualsFile = [subString '_task-' whichTask '_eegresidual' stretchMethod];
        save(fullfile(residualsFolder,residualsFile),'rEEGWithTS','rEEGWithFixed','beta','toExclude');
    end
    
    if ~doReg
        if doLaplacian
            % residualsFolder = ['./data/derivatives/eegresiduallap/sub-' participants{iParticipant}];
            residualsFolder = fullfile(dataFolder,'derivatives/eegresiduallap',subString);
        else
            % residualsFolder = ['./data/derivatives/eegresidual/sub-' participants{iParticipant}];
            residualsFolder = fullfile(dataFolder,'derivatives/eegresidual',subString);
        end
        if ~exist(residualsFolder,'dir'), mkdir(residualsFolder); end
        residualsFile = [subString '_task-' whichTask '_eegresidual_noreg'];
        save(fullfile(residualsFolder,residualsFile),'rEEGWithTS','rEEGWithFixed','beta','toExclude');
    end
end

% Save regression results
if doReg
    save(fullfile(resultsFolder,['rerp_results_' whichTask '_task_' stretchMethod '_' num2str(doLaplacian) '.mat']));
elseif ~doReg
    save(fullfile(resultsFolder, ['rerp_results_' whichTask '_task_noreg_' num2str(doLaplacian) '.mat']));
end
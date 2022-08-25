% Analysis in line with rest of temporal scaling paper

% Load preprocessed EEG and make ERPs
%
% Project: Temporal Scaling
% Other m-files required: EEGLAB, make_erp.m, find_artifacts.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% March 2020; Last revision: 02-Dec-2020

% Set raw and results folder (change as needed)
rawFolder = '/Users/chassall/Raw'; % From https://datadryad.org/stash/dataset/doi:10.5061/dryad.5vb8h
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

participants = {'02','03','04','05','06','07','08','10','11','12','13','14','15','16','18','19','20','21'};
nParticipants = length(participants);
judgeTimes = [0.512 0.64 0.8 1 1.25; 1.056 1.32 1.65 2.0625 2.5781; 1.6 2 2.5 3.125 3.9063];
subconditions = {'fast','slow'};
srate = 200;
doLaplacian = 0;
doReref = 0;

if doReref
   nElectrodes = 62; 
else
   nElectrodes = 64; 
end

%% Variables

% ERP windows
stimInterval = [-0.2,0.6];
respInterval = [-0.6,0.2];

% ERP baselines
stimBaseline = [-200,0];
respBaseline = [];

% ERP triggers
fastStimTriggers = {3};
slowStimTriggers = {4};
fastRespTriggers = {5};
slowRespTriggers = {6};

% ERP data
fastStimERP = nan(nParticipants,nElectrodes, round(srate * (stimInterval(2) - stimInterval(1))));
slowStimERP = nan(nParticipants,nElectrodes, round(srate * (stimInterval(2) - stimInterval(1))));
fastRespERP = nan(nParticipants,nElectrodes, round(srate * (respInterval(2) - respInterval(1))));
slowRespERP = nan(nParticipants,nElectrodes, round(srate * (respInterval(2) - respInterval(1))));

% Artifact proportions for each ERP of interest (see header for order)
allArtifacts = nan(nParticipants, 4);

meanRTs = [];

% Loop through participants
for iParticipant = 1:nParticipants
    
    % Load preprocessed data
    load(['sub-' participants{iParticipant}]);

    if doReref
        EEG = pop_reref(EEG,{'TP7','TP8'});
    end
    
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
    
    % Make ERPs
    [fastStimERP(iParticipant,:,:), allArtifacts(iParticipant,1)] = make_erp(EEG,fastStimTriggers,stimInterval,stimBaseline);
    [slowStimERP(iParticipant,:,:), allArtifacts(iParticipant,2)] = make_erp(EEG,slowStimTriggers,stimInterval,stimBaseline);
    [fastRespERP(iParticipant,:,:), allArtifacts(iParticipant,3)] = make_erp(EEG,fastRespTriggers,respInterval,respBaseline);
    [slowRespERP(iParticipant,:,:), allArtifacts(iParticipant,4)] = make_erp(EEG,slowRespTriggers,respInterval,respBaseline);
end

% Grand averages
fastStimGrandAve = squeeze(nanmean(fastStimERP));
slowStimGrandAve = squeeze(nanmean(slowStimERP));
fastRespGrandAve = squeeze(nanmean(fastRespERP));
slowRespGrandAve = squeeze(nanmean(slowRespERP));

% Save everything
if doReref
    save(fullfile(resultsFolder,'erp_results_dm_task_reref.mat'));
else
    save(fullfile(resultsFolder,'erp_results_dm_task.mat'));
end

%% rERP part
close all; clear all; rng(41); % for reproducibility

% Set raw and results folder (change as needed)
rawFolder = '/Users/chassall/Raw'; % From https://datadryad.org/stash/dataset/doi:10.5061/dryad.5vb8h
dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task4_dm/data';
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

participants = {'02','03','04','05','06','07','08','10','11','12','13','14','15','16','18','19','20','21'};
nParticipants = length(participants);
judgeTimes = [0.512 0.64 0.8 1 1.25; 1.056 1.32 1.65 2.0625 2.5781; 1.6 2 2.5 3.125 3.9063];
subconditions = {'fast','slow'};
srate = 200;
doLaplacian = 0;
doReref = 1; % Re-reference to the average of the mastoids?
if doReref
   nElectrodes = 62; 
else
    nElectrodes = 64;
end

% Method for stretching/compressing the scaled regressors
stretchMethod = 'box'; % 'box', 'triangle', 'cubic', or 'nearest'

% Regularization parameters
loadPrevReg = 1; % Re-use previous regularization result (saves time)
k = 10; % Number of fold for cross-validation
regChannels = {'Pz'}; % Do the cross-validation at these channels only
tsBreakpoints = [200 400]; % Since the betas will be appended, we need to know where the boundaries are
regtype = 'onediff'; % Use a finite difference approximation of the first derivative
lambdas = [0.001 0.01 0 1 10 100 1000 10000 100000 1000000 10000000]; % Reg. parameters to try
allCVErrors = nan(length(participants),nElectrodes,length(lambdas)); % To store cross-validation errors

% Stim-locked component
stimTimelimits = [-0.2,0.8]; % Time window
stimScalePntLimits = srate * stimTimelimits;
stimScalePnts = srate * stimTimelimits(1):srate *stimTimelimits(2);
numStimPnts = length(stimScalePnts);
stimBL = [-0.2 0];

% Response-locked component
respTimelimits = [-0.8,0.2]; % Time window
respScalePntLimits = srate * respTimelimits;
respScalePnts = srate * respTimelimits(1): srate *respTimelimits(2);
numRespPnts = length(respScalePnts);
respBL = [-0.8 -0.750]; % Baseline, in seconds

% Time-scaled component
scaleLength = 153; % In data points

% Load previous regularization results? (saves having to do it again)
if loadPrevReg && doReref
    load(fullfile(resultsFolder,'rerp_results_dm_task_reref_0.mat'));
elseif loadPrevReg
    load(fullfile(resultsFolder,'rerp_results_dm_task_0.mat'));
end

tsBreakpoints = [numStimPnts numStimPnts+numRespPnts];

%% Output variables

% Betas
allBetasStim = nan(length(participants),nElectrodes,numStimPnts);
allBetasResp = nan(length(participants),nElectrodes,numRespPnts);
allBetasTS  = nan(length(participants),nElectrodes,scaleLength);

% Error
allRSquared = nan(length(participants),nElectrodes);
allRSquaredFixed = nan(length(participants),nElectrodes);

% Residuals
allResiduals = [];
allStimRs = [];
allRespRs = [];

% Artifacts
allContinuousEEGArtifacts = nan(length(participants),1);
allContinuousBehaviouralArtifacts = nan(length(participants),1);

% Mean proportion of artifacts of either type
allContinuousArtifacts = nan(length(participants),1);

% Store design matrix
allX = {};

%% Participant Loop
% Loop through all participants, run the GLM, save the results
for iParticipant = 1:nParticipants
   
    % Load preprocessed data
    disp(participants{iParticipant});
    load(['sub-' participants{iParticipant}]);
    
    if doReref
    EEG = pop_reref(EEG,{'TP7','TP8'});
    end
    iRegChannels = eeg_chaninds(EEG,regChannels);
    
    % Flag EEG artifacts
    toExclude = [];
    isEEGArtifact = false(1,size(EEG.data,2));
    for a = 1:size(winrej,1)
        toExclude = [toExclude winrej(a,1):winrej(a,2)];
    end
    isEEGArtifact(toExclude) = true;
    
    % Make a copy of EEG that will eventually contain the continuous residual
    rEEGWithTS = EEG;
    rEEGWithFixed = EEG;
    
    stimX = sparse(zeros(size(EEG.data,2),numStimPnts));
    respX = sparse(zeros(size(EEG.data,2),numRespPnts));
    scalX = sparse(zeros(size(EEG.data,2),scaleLength));
    beta = nan(numStimPnts + numRespPnts +scaleLength,size(EEG.data,1));
    stimCount = 1;
    respCount = 1;
    for i = 1:length(EEG.event)
        thisLatency = round(EEG.event(i).latency);
        if EEG.event(i).type == 1
            stimX(thisLatency + stimScalePnts,:) = eye(numStimPnts);
            
            % Scaled component
            whichPnts = thisLatency:round(EEG.event(i+1).latency);
            thisBasis = imresize(eye(scaleLength),[length(whichPnts),scaleLength],stretchMethod);
            scalX(whichPnts,:) = thisBasis;
        elseif EEG.event(i).type == 2
            respX(thisLatency + respScalePnts,:) = eye(numRespPnts);
        end
    end
    % Truncate in case the regressors went past the end of the recording
    stimX = stimX(1:size(EEG.data,2),:);
    respX = respX(1:size(EEG.data,2),:);
    scalX = scalX(1:size(EEG.data,2),:);
    
    % Construct design matrix
    XNoTS = [stimX respX]; % Fixed-time components only (no time-scaling)
    XNoFixed = scalX;
    X = [stimX respX scalX]; % Fixed and scaled components
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
    %isTrialBehavArtifact = isBehaviouralArtifact & ~isEmptyRow;
    %allContinuousBehaviouralArtifacts(iParticipant) = mean(double(isTrialBehavArtifact));
    
    % isEitherArtifact = isTrialEEGArtifact | isTrialBehavArtifact;
    isEitherArtifact = isTrialEEGArtifact;
    isEitherArtifact(isEmptyRow) = [];
    allContinuousArtifacts(iParticipant) = mean(isEitherArtifact);
    
    
    isZero = sum(abs(X),2) == 0;
    lsmriterations = 400;
    iBad = [];
    for i = 1:size(winrej,1)
        iBad =  [iBad winrej(i,1):winrej(i,2)];
    end
    isBad = false(1,size(EEG.data,2));
    isBad(iBad) = true;
    %     for e = 1:size(EEG.data,1)
    %         thisData = EEG.data(e,:);
    %         thisX = X;
    %         thisData(isZero' | isBad) = [];
    %         thisX(isZero | isBad',:) = [];
    %         [beta(:,e),ISTOP,ITN] = lsmr(thisX,double(thisData'),[],10^-8,10^-8,[],lsmriterations);
    %     end
    
    
    
    toExclude = isZero' | isBad;
    X(toExclude,:) = [];
    data(:,toExclude) = [];
    XFixed(toExclude,:) = [];
    hasFixed(toExclude) = [];

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
    plot(allErrors(iRegChannels,:)); drawnow();
    [M,I] = min(allErrors,[],2);
    bestLambda = lambdas(I(iRegChannels));
    % bestLambda = 1000; % Hard code for now
    % regtype = 'onediff';
    
    % Do regression
    pDM = pinv_reg(X,bestLambda,regtype,tsBreakpoints);
    pDMFixed = pinv_reg(XFixed,bestLambda,regtype,tsBreakpoints(1));
    data = double(data);
    for c = 1:EEG.nbchan
        beta(:,c) = pDM * data(c,:)';
        betaFixed(:,c) = pDMFixed * data(c,:)';
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
        
    % Compute residuals
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
    if doLaplacian
        residualsFolder = fullfile(dataFolder,'derivatives','eegresiduallap',['sub-' participants{iParticipant}]);
    else
        residualsFolder = fullfile(dataFolder,'derivatives','eegresidual',['sub-' participants{iParticipant}]);
    end
    if ~exist(residualsFolder,'dir'), mkdir(residualsFolder); end
    if doReref
        residualsFile = ['sub-' participants{iParticipant} '_task-dm_eegresidual_reref'];
    else
        residualsFile = ['sub-' participants{iParticipant} '_task-dm_eegresidual'];
    end
    save(fullfile(residualsFolder,residualsFile),'rEEGWithTS','rEEGWithFixed','beta','toExclude');
end

% Save regression results
if doReref
    save(fullfile(resultsFolder,['rerp_results_dm_task_reref_' num2str(doLaplacian) '.mat']));
else
    save(fullfile(resultsFolder,['rerp_results_dm_task_' num2str(doLaplacian) '.mat']));
end

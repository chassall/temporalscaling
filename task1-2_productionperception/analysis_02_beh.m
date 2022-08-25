% Load and analyze behavioural data
%
% Project: Temporal Scaling

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% March 2020; Last revision: 11-Jul-2022

participants = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};
nParticipants = length(participants);

% Data and results folder. Change as needed.
dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/task1-2_productionperception/data';
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

% Output Variables
allDoWindows = nan(nParticipants,3);
allDoIntervals = nan(nParticipants,3);
allDoAdjustments = nan(nParticipants,3,3); % p X condition  X FB (early, on time, late)
allJudgeYesProb = nan(nParticipants,5,3); % p X subcondition X condition

doLim = [0.2 0.8*5; 0.2 1.65*5; 0.2 2.5*5]; % lower/upper limits to cut out really bad trials

for iParticipant = 1:nParticipants
    
    %% Load behavioural data
    % Data columns are as follows (see behavioural README for more detail):
    % block, trial, pre-cue jitter, RT (do task), computer RT (judge task),
    % probe condition (very early, early, etc.), response mapping (judge
    % task), response, judgement (judge task), second RT (judge task),
    % pre-feedback jitter, feedback condition, point total, response margin
    % (do task)
    subString = ['sub-' participants{iParticipant}];
    % behFolder = ['./data/rawdata/sub-' participants{iParticipant} '/beh'];
    behFolder = fullfile(dataFolder,subString,'beh');
    behFile = [subString '_task-temporalscaling_beh.tsv'];
    
    % Load tab-separated behavioural data
    opts = delimitedTextImportOptions("NumVariables", 14);
    opts.Delimiter = "\t";
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    theseData = readtable(fullfile(behFolder,behFile), opts);
    theseData = table2array(theseData); % Do this because we originally worked with an array

    % Feedback conditions
    % 1 incorrect (judge)
    % 2 correct (judge)
    % 3 do early (do)
    % 4 do correct (do)
    % 5 do late (do)
    
    %% Isolate data for each task (do, judge)
    doTrials = ismember(theseData(:,1),[1 2 3]);
    doData  = theseData(doTrials,:);
    judgeTrials = ismember(theseData(:,1),[4 5 6]);
    judgeData = theseData(judgeTrials,:);
    
    %% Do the "judge" analysis (Temporal Prediction)
    
    judgeSubconditions = [1 2 3 4 5]; % v.early, early, on time, late, v. late
    probYesResponse = [];
    judgeBlocks = [4 5 6];
    
    % Loop through each condition (short, medium, long)
    % and subcondition (v. early, early, on time, late, v. late).
    % Compute the proportion of "yes" responses.
    for iCondition = 1:3
        isCondition = (judgeData(:,1) == judgeBlocks(iCondition));
        for iSubcondition = 1:length(judgeSubconditions)
            isSubcondition = judgeSubconditions(iSubcondition);
            iTrial = isCondition & judgeData(:,6) == isSubcondition;
            sumYes = sum(judgeData(iTrial,9) == 1); % Yes
            sumNo = sum(judgeData(iTrial,9) == 0); % No
            probYesResponse(iSubcondition, iCondition) = 100 * sumYes / (sumYes + sumNo);
        end
        
    end
    allJudgeYesProb(iParticipant,:,:) = probYesResponse;
   
    %% "Do" task (Temporal Production)
    
    doBlocks = [1 2 3]; % short, medium, long
    doFeedbackConditions = [3 4 5]; % early, correct, late
   

    % Loop through conditions and compute mean window size and RT
    for iCondition = 1:3
        

        isCondition = (doData(:,1) == iCondition);
        isBad = doData(:,4) < doLim(iCondition,1) | doData(:,4) > doLim(iCondition,2);
        allDoWindows(iParticipant,iCondition) = mean(1000*doData(isCondition,14));
        allDoIntervals(iParticipant,iCondition) = mean(1000*doData(isCondition & ~isBad,4));
    end
    
    % Loop through feedback types (early, correct, late)
    for fbType = 1:length(doFeedbackConditions)
        % Loop through conditions (short, medium, long)
        for iCondition = 1:length(doBlocks)
           isCondition =  doData(:,1) == doBlocks(iCondition);
           % Isolate condition data
           thisData = doData(isCondition,:); 
           % Get condition RTs
           theseRTs = thisData(:,4);
           % Get "trial n" feedback
           theseFeedbacks = thisData(1:end-1,12);
           % Get "trial n+1" RT adjustment
           theseAdjustments = theseRTs(2:end) - mean(theseRTs);
           % Get "trial n+1" adjustment for this "trial n" feedback
           isThisFeedback = theseFeedbacks == doFeedbackConditions(fbType);
           allDoAdjustments(iParticipant,iCondition,fbType) = nanmean(theseAdjustments(isThisFeedback));
        end
    end  
end % End participant loop

% Save results
save(fullfile(resultsFolder,'beh_results_production_and_perception_tasks.mat'));


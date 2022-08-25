% Behavioural analysis for the Prediction task
%

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% 25-Aug-2022

close all; clear all; clc;
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

%% Prediction Task

% Column number	Variable
% 1	Subject number
% 2	Block number
% 3	Block type (0 = Non-Informative, 1 = Rhythmic, 2 = Repeated-Interval, 3 = Random)
% 4	Trial number
% 5	Separator (-1)
% 6	Trial type (same as column 3)
% 7	Stream SOA (1 = short, 2 = long)
% 8	Target SOA (1 = short, 2 = long, 0 = catch)
% 9	Cue validity (1 = valid, 0 = invalid, 2 = catch)
% 10	Number of stream stimuli
% 11	Response
% 12	RT
% 13	Jitter of target from mean of SOA distribution (+-4 steps of 50 msec)

ps = [1:12 14 16:21]; % 13, 15 missing on datadryad
participants = num2str(ps,'%0.2d,');
participants(end) = [];
participants = strsplit(participants,',');

behFolder = 'E:\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_BIDS\task3_prediction\data\beh';

allMeans = [];
for i = 1:length(participants)
   thisFile = fullfile(behFolder,['entrainRhythm1B_EEG_subject1' participants{i} '.txt']);
   thisData = load(thisFile); 
   
   isRhythmic = thisData(:,3) == 1;
   isRepeated = thisData(:,3) == 2;
   isValid = thisData(:,9) == 1;
   isShort = thisData(:,8) == 1;
   isLong = thisData(:,8) == 2;
   
   cond1 = isRhythmic & isValid & isShort;
   cond2 = isRhythmic & isValid & isLong;
   cond3 = isRepeated & isValid & isShort;
   cond4 = isRepeated & isValid & isLong;
   
   allMeans(i,1) = mean(thisData(cond1,12));
   allMeans(i,2) = mean(thisData(cond2,12));
   allMeans(i,3) = mean(thisData(cond3,12));
   allMeans(i,4) = mean(thisData(cond4,12));
end

save(fullfile(resultsFolder,'beh_prediction.mat'));
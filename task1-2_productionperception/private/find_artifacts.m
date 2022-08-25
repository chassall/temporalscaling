%FIND_ARTIFACTS Flags bag epochs
%   Takes an EEGLAB data structure that is epoched and flags bad epochs
%   based on the sample-to-sample change, or on the overall change within
%   the epoch.
%
% Syntax:  [isArtifact, isArtifactsCT] = find_artifacts(EEG,threshold)
%
% Inputs:
%    EEG - EEGLAB data structure where the EEG has been epoched
%    threshold (optional) - max sample-to-sample diff and max diff in epoch
%
% Outputs:
%    isArtifact - Vector of artifact flags (1 = artifact, 0 = OK)
%    isArtifactsCT - Channels X Trials artifacts (1 = artifact, 0 = OK)
%
% Example: 
%    badTrialIndex = findArtifacts(icaEEG,1000); % 1000 uV cutoff
%    [badTrialIndex, badTrialsByChannel] = findArtifacts(icaEEG,[40,125]);
%    badTrialIndex = findArtifacts(icaEEG);
%
% Author: Cameron Hassall, Department of psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% September 2019; Last revision: 03-Jun-2020

function [isArtifact, isArtifactsCT] = find_artifacts(EEG, maxMin, level, step, lowest)

if nargin == 1
    maxMin = 150;
    level = 150;
    step = 40;
    lowest = 0.1;
end

if length(nargin) == 2
    isArtifactsCT = abs((max(EEG.data,[],2) -  min(EEG.data,[],2))) > maxMin;
else
    isMaxMinArtifact = abs((max(EEG.data,[],2) -  min(EEG.data,[],2))) > maxMin;
    isLevelArtifact = max(abs(EEG.data),[],2) > level;
    isStepArtifact = any(diff(EEG.data,[],2) > step,2);
    isLowestArtifact = all(abs(EEG.data) < lowest,2);
    isArtifactsCT = isMaxMinArtifact | isLevelArtifact | isStepArtifact | isLowestArtifact;
    isArtifactsCT = squeeze(isArtifactsCT);
end
isArtifact = logical(squeeze(any(isArtifactsCT,1)));

end

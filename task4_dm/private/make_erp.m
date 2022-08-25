function [erp,artifactProportion, times, data, isArtifactsCT] = make_erp(EEG,triggers,interval,baseline,artifactSettings)
%MAKE_ERP Makes an ERP

%maxMin = 100;
%level = 100;
%step = 40;
%lowest = 0.1;
%     maxMin = 400;
%     level = 400;
%     step = 40;
%     lowest = 0.1;

if nargin < 4
    baseline = [];
    maxMin = 100;
    level = 100;
    step = 40;
    lowest = 0.1;
elseif nargin < 5
    maxMin = 100;
    level = 100;
    step = 40;
    lowest = 0.1;
else
    maxMin = artifactSettings.maxMin;
    level = artifactSettings.level;
    step = artifactSettings.step;
    lowest = artifactSettings.lowest;
end

%   Detailed explanation goes here
% Make epochs
EEG = pop_epoch(EEG,triggers,interval);

% Baseline correction
if ~isempty(baseline)
    EEG = pop_rmbase(EEG,baseline);
end

% Artifact rejection
[isArtifact, isArtifactsCT] = find_artifacts(EEG, maxMin, level, step, lowest);
artifactProportion = mean(isArtifactsCT,2);

% Make ERP
erp = mean(EEG.data(:,:,~isArtifact),3);

times = EEG.times;

data = EEG.data;
end


% Supplementary Figure 9, Variance Inflation Factor
%
% Other m-files required: subtightplot, plot_settings

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% 25-Aug-2022

% Other files neededed: unfold toolbox

init_unfold(); % Needed for uf_vif;
close all; clear all; clc;

resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';
tiffFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/figures/tiffs';
tiffRes = '-r600';
whichTask = 'dm';

switch whichTask
    
    case 'production'
        load(fullfile(resultsFolder,'rerp_results_production_task_box_0.mat'));
        stimTimes = stimTimelimits(1):1/srate:stimTimelimits(2);
        respTimes = respTimelimits(1):1/srate:respTimelimits(2);
        scaleTimes = linspace(1,100,scaleLength);
        titles = {'Cue','Response','Scaled'};
    case 'perception'
        load(fullfile(resultsFolder,'rerp_results_perception_task_box_0.mat'));
        stimTimes = stimTimelimits(1):1/srate:stimTimelimits(2);
        respTimes = respTimelimits(1):1/srate:respTimelimits(2);
        scaleTimes = linspace(1,100,scaleLength);
        titles = {'Cue','Probe','Scaled'};
    case 'prediction'
        load(fullfile(resultsFolder,'rerp_and_erp_results_prediction_task_0.mat'));
        stimTimes = warningTimes;
        respTimes = targetTimes;
        scaleTimes = linspace(1,100,scaleLength);
        titles = {'Cue','Target','Scaled'};
    case 'dm'
        load(fullfile(resultsFolder,'rerp_results_dm_task_reref_0.mat'));
        stimTimes = stimTimelimits(1):1/srate:stimTimelimits(2);
        respTimes = respTimelimits(1):1/srate:respTimelimits(2);
        scaleTimes = linspace(1,100,scaleLength);
        titles = {'Cue','Decision','Scaled'};
end

% load('../results/rerp_results_perception_task_box_0.mat');
% load('../results/rerp_results_production_task_box_0.mat');
% stimTimes = stimTimelimits(1):1/srate:stimTimelimits(2);
% respTimes = respTimelimits(1):1/srate:respTimelimits(2);
% scaleTimes = linspace(1,100,scaleLength);

%%
allVIF = [];
for p = 1:length(allX)
   thisX = allX{p};
   thisVIF = uf_vif(thisX);
   allVIF(p,:) = thisVIF;
end

meanVIF = mean(allVIF,1);
stdVIF = std(allVIF,[],1);

%% Plot
plot_settings;
stimVIF = meanVIF(1:tsBreakpoints(1));
respVIF = meanVIF((tsBreakpoints(1)+1):tsBreakpoints(2));
scaledVIF = meanVIF((tsBreakpoints(2)+1):end);

stimCI = stdVIF(1:tsBreakpoints(1));
respCI = stdVIF((tsBreakpoints(1)+1):tsBreakpoints(2));
scaledCI = stdVIF((tsBreakpoints(2)+1):end);

axs = [];
makefigure(19,4);
% axs(1) = subplot(1,3,1); plot(stimTimes,stimVIF,'k-','LineWidth',plotLineWidth); box off; xlabel('Time (ms)'); ylabel('VIF');
% axs(2) = subplot(1,3,2); plot(respTimes,respVIF,'LineWidth',plotLineWidth); box off; xlabel('Time (ms)'); ylabel('VIF');
% axs(3) = subplot(1,3,3); plot(scaleTimes,scaledVIF,'LineWidth',plotLineWidth); box off; xlabel('Interval proportion (%)'); ylabel('VIF');

axs(1) = subplot(1,3,1); boundedline(stimTimes,stimVIF,stimCI); xlabel('Time (ms)'); ylabel('VIF'); title(titles{1});
axs(2) = subplot(1,3,2); boundedline(respTimes,respVIF,respCI); xlabel('Time (ms)'); ylabel('VIF'); title(titles{2});
axs(3) = subplot(1,3,3); boundedline(scaleTimes,scaledVIF,scaledCI);  xlabel('Interval proportion (%)'); ylabel('VIF'); title(titles{3});

linkaxes(axs,'y');

switch whichTask
    
    case 'prediction'
        ylim([0 100]);
    case 'dm'
        ylim([0 50]);
end

print(fullfile(tiffFolder,['figure_s08_vif_' whichTask '.tiff']),'-dtiff',tiffRes);

%% Plot
stimVIF = allVIF(:,1:tsBreakpoints(1));
respVIF = allVIF(:,(tsBreakpoints(1)+1):tsBreakpoints(2));
scaledVIF = allVIF(:,(tsBreakpoints(2)+1):end);

axs = [];
makefigure();
axs(1) = subplot(1,3,1); plot(stimTimes,stimVIF');
axs(2) = subplot(1,3,2); plot(respTimes,respVIF');
axs(3) = subplot(1,3,3); plot(scaleTimes,scaledVIF');

linkaxes(axs,'y');
% Supplementary Figure 4, ERPssubtightplot
%
% Other m-files required: subtightplot, plot_settings

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% 25-Aug-2022

%% Prodiction and Perception tasks
close all; clear all;
resultsFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\results';
tiffRes = '-r600';
load(fullfile(resultsFolder,'erp_results_production_and_perception_tasks.mat'));
outputFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\figures\tiffs';

channel = 'FCz';
iChannel = eeg_chaninds(EEG,{channel});

plot_settings;
% makefigure(33.87, 16);
makefigure(18.26,4);
subplot = @(m,n,p) subtightplot (m, n, p, [0.2 0.08], [0.2 0.1], [0.05 0.02]);

ax(1) = subplot(1,2,1);
plot(doShortBeepInterval(1):(1/srate):doShortBeepInterval(2)-(1/srate),doShortBeepGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(1,:)); hold on;
plot(doMediumBeepInterval(1):(1/srate):doMediumBeepInterval(2)-(1/srate),doMediumBeepGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(2,:));
plot(doLongBeepInterval(1):(1/srate):doLongBeepInterval(2)-(1/srate),doLongBeepGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(5,:));
title('Cue'); legend('Short','Medium','Long', 'box','off'); xlim([-0.2,2.5]);

ax(2) = subplot(1,2,2);
plot(doShortResponseInterval(1):(1/srate):doShortResponseInterval(2)-(1/srate),doShortResponseGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(1,:)); hold on;
plot(doMediumResponseInterval(1):(1/srate):doMediumResponseInterval(2)-(1/srate),doMediumResponseGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(2,:));
plot(doLongResponseInterval(1):(1/srate):doLongResponseInterval(2)-(1/srate),doLongResponseGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(5,:));
title('Response'); legend('Short','Medium','Long', 'box','off','Location','southwest'); xlim([-2.5,0.2]); 

for i = 1:2
    ax(i).XLabel.String = 'Time (s)';
    ax(i).YLabel.String = 'Voltage (\muV)';
    ax(i).FontSize = plotFontSize;
    ax(i).FontName = plotFontName;
    ax(i).Box = 'off';
    t = text(ax(i).XLim(1)+0.02*diff(ax(i).XLim),ax(i).YLim(2),channel,'Parent', ax(i));
    t.FontSize = plotFontSize;
end

print(fullfile(outputFolder,'sfig_04_production_erps.tiff'),'-dtiff',tiffRes);

makefigure(18.26,4);
ax(3) = subplot(1,2,1);
plot(judgeShortBeepInterval(1):(1/srate):judgeShortBeepInterval(2)-(1/srate),judgeShortBeepGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(1,:)); hold on;
plot(judgeMediumBeepInterval(1):(1/srate):judgeMediumBeepInterval(2)-(1/srate),judgeMediumBeepGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(2,:));
plot(judgeLongBeepInterval(1):(1/srate):judgeLongBeepInterval(2)-(1/srate),judgeLongBeepGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(5,:));
title('Cue'); legend('Short','Medium','Long', 'box','off'); xlim([-0.2,2.5]);

ax(4) = subplot(1,2,2);
plot(judgeShortProbeInterval(1):(1/srate):judgeShortProbeInterval(2)-(1/srate),judgeShortProbeGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(1,:)); hold on;
plot(judgeMediumProbeInterval(1):(1/srate):judgeMediumProbeInterval(2)-(1/srate),judgeMediumProbeGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(2,:));
plot(judgeLongProbeInterval(1):(1/srate):judgeLongProbeInterval(2)-(1/srate),judgeLongProbeGrandAve(iChannel,:),'LineWidth',plotLineWidth,'Color',eventColours(5,:));
title('Probe'); legend('Short','Medium','Long', 'box','off','Location','southwest'); xlim([-2.5,0.2]);

for i = 2:4
    ax(i).XLabel.String = 'Time (s)';
    ax(i).YLabel.String = 'Voltage (\muV)';
    ax(i).FontSize = plotFontSize;
    ax(i).FontName = plotFontName;
    ax(i).Box = 'off';
    t = text(ax(i).XLim(1)+0.02*diff(ax(i).XLim),ax(i).YLim(2),channel,'Parent', ax(i));
    t.FontSize = plotFontSize;
end

% print([outputFolder 'figure_s01_perception_erps.tiff'],'-dtiff','-r300');
print(fullfile(outputFolder,'sfig_04_perception_erps.tiff'),'-dtiff',tiffRes);

%% Prediction Task

close all; clear all;
resultsFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\results';
outputFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\figures\tiffs';
tiffRes = '-r600';
load(fullfile(resultsFolder,'rerp_and_erp_results_prediction_task_0.mat'));

plot_settings;
makefigure(18.26,4);
channel = 'FCz';
c = eeg_chaninds(EEG,{channel});

subplot = @(m,n,p) subtightplot (m, n, p, [0.2 0.08], [0.18 0.1], [0.05 0.02]);

ax(1) = subplot(1,2,1);
plot(wEEG.times./1000,squeeze(mean(erps(:,1,c,:),1)),'LineWidth',plotLineWidth,'Color',eventColours(1,:));
hold on
plot(longCueEEG.times./1000,squeeze(mean(longCueERPs(:,c,:),1)),'LineWidth',plotLineWidth,'Color',eventColours(5,:));
l = legend('Short','Long');
l.Box = 'off';
l.Location = 'SouthWest';
title('Cue');
ylabel('Voltage (\muV)');

ax(2) = subplot(1,2,2);
plot(tEEG.times./1000,squeeze(mean(erps(:,3,c,:),1)),'LineWidth',plotLineWidth,'Color',eventColours(1,:));
hold on;
plot(longTargetEEG.times./1000,squeeze(mean(longTargetERPs(:,c,:),1)),'LineWidth',plotLineWidth,'Color',eventColours(5,:));
l = legend('Short','Long');
l.Box = 'off';
l.Location = 'SouthWest';
title('Target');
ylabel('Voltage (\muV)');

for i = 1:2
    ax(i).XLabel.String = 'Time (s)';
    ax(i).YLabel.String = 'Voltage (\muV)';
    ax(i).FontSize = plotFontSize;
    ax(i).FontName = plotFontName;
    ax(i).Box = 'off';
    t = text(ax(i).XLim(1)+0.02*diff(ax(i).XLim),ax(i).YLim(2),channel,'Parent', ax(i));
    t.FontSize = plotFontSize;
end

ax(1).XLim = [-0.2,1.3];
ax(2).XLim = [-1.3,0.2];

print(fullfile(outputFolder,'sfig_04_prediction_erps.tiff'),'-dtiff',tiffRes);

%% DM Task

close all; clear all;
resultsFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\results';
outputFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\figures\tiffs';
tiffRes = '-r600';
load(fullfile(resultsFolder,'erp_results_dm_task.mat'));

plot_settings;
makefigure(18.26,4);
channel = 'Pz';
stimTimes = stimInterval(1):1/EEG.srate:stimInterval(2)-1/EEG.srate;
respTimes = respInterval(1):1/EEG.srate:respInterval(2)-1/EEG.srate;
c = eeg_chaninds(EEG,{channel});

subplot = @(m,n,p) subtightplot (m, n, p, [0.2 0.08], [0.18 0.1], [0.05 0.02]);

ax(1) = subplot(1,2,1);
plot(stimTimes,squeeze(mean(fastStimERP(:,c,:),1)),'LineWidth',plotLineWidth,'Color',eventColours(1,:));
hold on
plot(stimTimes,squeeze(mean(slowStimERP(:,c,:),1)),'LineWidth',plotLineWidth,'Color',eventColours(5,:));
l = legend('Fast','Slow');
l.Box = 'off';
l.Location = 'NorthWest';
title('Cue');
ylabel('Voltage (\muV)');

ax(2) = subplot(1,2,2);
plot(respTimes,squeeze(mean(fastRespERP(:,c,:),1)),'LineWidth',plotLineWidth,'Color',eventColours(1,:));
hold on;
plot(respTimes,squeeze(mean(slowRespERP(:,c,:),1)),'LineWidth',plotLineWidth,'Color',eventColours(5,:));
l = legend('Fast','Slow');
l.Box = 'off';
l.Location = 'Best';
title('Decision');
ylabel('Voltage (\muV)');

for i = 1:2
    ax(i).XLabel.String = 'Time (s)';
    ax(i).YLabel.String = 'Voltage (\muV)';
    ax(i).FontSize = plotFontSize;
    ax(i).FontName = plotFontName;
    ax(i).Box = 'off';
    t = text(ax(i).XLim(1)+0.02*diff(ax(i).XLim),ax(i).YLim(2),channel,'Parent', ax(i));
    t.FontSize = plotFontSize;
end

print(fullfile(outputFolder,'sfig_04_dm_erps.tiff'),'-dtiff',tiffRes);

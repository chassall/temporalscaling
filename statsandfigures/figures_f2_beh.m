% Loads and plots the behavioural data
%
% Other m-files required: 

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% Dec 2020; Last revision: 22-Dec-2021

close all; clear all; clc;
allInOne = 0; plotHeight = 3;

resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';
tiffFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/figures/tiffs';
tiffRes = '-r600';

%% Load production/perception behavioural results
load(fullfile(resultsFolder,'beh_results_production_and_perception_tasks.mat'));
clear subplot;
plot_settings;

% subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.14], [0.2 0.1], [0.1 0.1]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

if allInOne
makefigure(18.3,15);
subplot(2,2,1);
else
makefigure(5.36,plotHeight);
end
[h2, stats2] = notBoxPlot(allDoIntervals,'interval','tInterval');

disp('"Do" task stats (RT)');
disp([stats2.mu]);
disp([stats2.mu] - [stats2.interval]);
disp([stats2.mu] + [stats2.interval]);

formatNBP(h2,plotColours);
hold on;
yline(800,':','Color',plotColours(1,:));
yline(1650,':','Color',plotColours(2,:));
yline(2500,':','Color',plotColours(3,:));
% ylim([0 3000]);

ax = gca;
ax.Box = 'off';
ax.FontSize = plotFontSize;
ax.FontName = plotFontName;
ax.XLabel.String = 'Trial type';
ax.YLabel.String = 'Response time (ms)';
ax.XTickLabel = {'Short','Medium','Long'};

if ~allInOne
print(fullfile(tiffFolder,'fig_2_productionbeh.tiff'),'-dtiff','-r600');
end

% "Judge" task
if allInOne
subplot(2,2,2);
else
    makefigure(5.36,plotHeight);
end
meanP = squeeze(mean(allJudgeYesProb,1));
stdP = squeeze(std(allJudgeYesProb,[],1));
tVal = abs(tinv(0.025,nParticipants-1));
ciP = tVal * stdP ./ sqrt(nParticipants);
logTimes = 10*log([1.25^-2 1.25^-1 1.25^0 1.25^1 1.25^2]);

disp('Perception Task');
disp(meanP);
disp(meanP-ciP);
disp(meanP+ciP);
allPs = [];
for i = 1:3
    allPs(i) = plot(logTimes,meanP(:,i),'Color',plotColours(i,:),'LineWidth',plotLineWidth);
    hold on;
    plot(logTimes,meanP(:,i),'o','MarkerFaceColor',plotColours(i,:),'MarkerSize',4);
    errorbar(logTimes,meanP(:,i),ciP(:,i),'Color',plotColours(i,:));
end

l = legend(allPs,'Short','Medium','Long');
l.Box = 'off';
l.Position = [0.5791    0.6250    0.4262    0.4044];
ylim([0 100]);
xlabel('Probe offset');
ylabel('"Yes" resp. (%)');
ax = gca;
ax.Box = 'off';
ax.FontSize = plotFontSize;
ax.FontName = plotFontName;
ax.XTick = logTimes;
ax.XTickLabel = {'Very early','Early','On time','Late','Very late'};
ax.XTickLabelRotation = -35;

if ~allInOne
% print('./tiffs/fig_2_perceptionbeh.tiff','-dtiff','-r600');
print(fullfile(tiffFolder,'fig_2_perceptionbeh.tiff'),'-dtiff','-r600');
end

%% Prediction Task
% load('../results/beh_prediction.mat');
load(fullfile(resultsFolder,'beh_prediction.mat'));

if allInOne
subplot(2,2,3);
else
makefigure(5.36,plotHeight);
end
[nbp, stats] = notBoxPlot(allMeans,'interval','tInterval');
formatNBP(nbp);
xticklabels({'Short','Long','Short','Long'});
ylabel('Reaction time (s)');
ylim([0.15 0.35]);
text(0.8,0.34,'Rhythmic','FontSize',8);
text(2.8,0.34,'Repeated','FontSize',8);

ax = gca;
ax.Box = 'off';
ax.FontSize = plotFontSize;
ax.FontName = plotFontName;

disp('Prediction Task');
[stats.mu]
[stats.mu] - [stats.interval]
[stats.mu] + [stats.interval]

if ~allInOne
% print('./tiffs/fig_2_predictionbeh.tiff','-dtiff','-r600');
print(fullfile(tiffFolder,'fig_2_predictionbeh.tiff'),'-dtiff','-r600');
end

%% DM Task
% load('../results/scalefactor_dm.mat','allRTs');
load(fullfile(resultsFolder,'scalefactor_dm.mat'));

if allInOne
subplot(2,2,4);
else
    makefigure(5.36,plotHeight);
end
histogram(allRTs);
xlabel('Response time (s)');
ylabel('Trial count');
xline(mean(allRTs),'r');

disp('Decision Making');
meanRT = mean(allRTs);
disp(meanRT)

text(0.36,550,'Fast','Color','r','FontSize',8);
text(0.87,550,'Slow','Color','r','FontSize',8);

ax = gca;
ax.Box = 'off';
ax.FontSize = plotFontSize;
ax.FontName = plotFontName;

if ~allInOne
% print('./tiffs/fig_2_dmbeh.tiff','-dtiff','-r600');
print(fullfile(tiffFolder,'fig_2_dmbeh.tiff'),'-dtiff','-r600');
end

%% Print

if allInOne
% print('./tiffs/sfig_2_beh.tiff','-dtiff',tiffRes);
print(fullfile(tiffFolder,'sfig_2_beh.tiff'),'-dtiff','-r600');
end
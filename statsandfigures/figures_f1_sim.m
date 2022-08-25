% Shows simulated components, draws a sample design matrix
%
% Other m-files required: subtightplot

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% June 2020; Last revision: 25-Aug-2022

close all; clear all; clc;

% Set results and output folders (change as needed)
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';
tiffFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/figures/tiffs';
tiffRes = '-r600';

% Plot settings
plotLineWidth = 1;
plotColours = cbrewer('qual','Dark2',5);
eventColours = cbrewer('qual','Set1',5);
eventColours = eventColours([1,2,5],:);
axisSettings.fontSize = 8;
axisSettings.fontName = 'Arial';

% Pick a simulation to plot
whichSimulation = 'trfparam'; % 'standard','cnv','flat','pulse','narrowpulse','boxcar','trf','trfparam','trfellapsedt'
load(fullfile(resultsFolder,['task1-2_simulationresults_' whichSimulation '.mat']));

%% Fig. 1a
gap = [0.08 0.01];
marg_h = [0.01 0.01];
marg_w = [0.01 001];
subplot = @(m,n,p) subtightplot (m, n, p,gap,marg_h,marg_w);
%   Input arguments (defaults exist):
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

thickPlotLineWidth = 1.5;
thinPlotLineWidth = 1; 
commonYLim = [-2 4];
commonWindow = 3000:4500;
makefigure(5.59,9.53);

subplot(6,1,1);
plot(0:400,stim_evoked,'LineWidth',thickPlotLineWidth,'Color',plotColours(1,:)); 
xlim([0 intervals(3)+100]);
set_axis(gca,axisSettings);
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

subplot(6,1,2);
% plot(-400:100,button_press,'LineWidth',thickPlotLineWidth,'Color',eventColours(2,:)); 
plot((intervals(1)-400):(intervals(1)+100),button_press,'LineWidth',thickPlotLineWidth,'Color',plotColours(2,:)); 
hold on;
plot((intervals(2)-400):(intervals(2)+100),button_press,'LineWidth',thickPlotLineWidth,'Color',plotColours(2,:)); 
plot((intervals(3)-400):(intervals(3)+100),button_press,'LineWidth',thickPlotLineWidth,'Color',plotColours(2,:)); 

xlim([0 intervals(3)+100]);
set_axis(gca,axisSettings);
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

subplot(6,1,3);
for i = 1:length(intervals)
    if strcmp(whichSimulation,'standard') 
    plot(0:intervals(i),tscs{i},'LineWidth',thickPlotLineWidth,'Color',plotColours(3,:)); hold on;
    else
    plot(0:intervals(i),tscs{i}+i,'LineWidth',thickPlotLineWidth,'Color',plotColours(3,:)); hold on;
    end
end
set_axis(gca,axisSettings);
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
% ax.YLim = commonYLim;

subplot(6,1,4);
plot(EEGnoise(commonWindow),'k','LineWidth',thinPlotLineWidth);
set_axis(gca,axisSettings);
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.YLim = commonYLim;
ax.XLim = [1 length(commonWindow)];

subplot(6,1,6);
plot(EEGdat(commonWindow),'k','LineWidth',thinPlotLineWidth);
set_axis(gca,axisSettings);
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.YLim = commonYLim;
ax.XLim = [1 length(commonWindow)];

print(fullfile(tiffFolder, ['fig_01a_simcomponents_' whichSimulation '.tiff']),'-dtiff',tiffRes);

%% Fig. 1b 
eegPnts = 1:30000;
makefigure(10,2);
plot(EEGdat(eegPnts),'k');
ax = gca;
ax.Box = 0;
ax.XAxis.Visible = 0;
ax.YAxis.Visible = 0;
plotLim = ax.YLim;
print(fullfile(tiffFolder,['fig_01b_simeeg_' whichSimulation '.tiff']),'-dtiff',tiffRes);

%% Fig. 1b
eegPnts = 1:30000;
makefigure(10,2);
plot(residuals(eegPnts),'k');
ax = gca;
ax.Box = 0;
ax.XAxis.Visible = 0;
ax.YAxis.Visible = 0;
ax.YLim = plotLim;
print(fullfile(tiffFolder,['fig_01b_simresiduals_' whichSimulation '.tiff']),'-dtiff',tiffRes);

%% Fig. 1b (unused)
gap = [0.01 0.01];
marg_h = [0.05 0.05];
marg_w = [0.01 0.01];
subplot = @(m,n,p) subtightplot (m, n, p,gap,marg_h,marg_w);
%   Input arguments (defaults exist):
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 


commonWindow = 2500:14000;
eegLim = [-5 10];
fig = makefigure(7.98,9.25);
testX = X(commonWindow,:);

subplot(1,4,1);
plot(EEGdat(commonWindow),1:length(commonWindow),'Color','k','LineWidth',plotLineWidth)
ax = gca;
ax.Box = 'off';
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.YLim = [1 length(commonWindow)];
ax.YDir = 'reverse';
ax.XLim = eegLim;

subplot(1,4,2);
cspy(X(commonWindow,:),'MarkerSize', 6,'Levels',1);
set_axis(gca,axisSettings);
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

subplot(1,4,4);
plot(EEGnoise(commonWindow),1:length(commonWindow),'Color','k','LineWidth',plotLineWidth)
ax = gca;
ax.Box = 'off';
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.YDir = 'reverse';
ax.YLim = [1 length(commonWindow)];
ax.XLim = eegLim;

pause(0.1);

%% Fig. 1c
gap = [0.16 0.12];
marg_h = [0.12 0.06];
marg_w = [0.1 0.06];
subplot = @(m,n,p) subtightplot (m, n, p,gap,marg_h, marg_w);
%   Input arguments (defaults exist):
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

thickPlotLineWidth = 1;
thinPlotLineWidth = 1; 
commonYLim = [-2 4];
commonWindow = 2800:4200;
makefigure(8,8);

subplot(3,1,1);
plot(0:0.001:0.999,beta(1:1000),'LineWidth',thickPlotLineWidth,'Color',plotColours(1,:)); 
set_axis(gca,axisSettings);
ax = gca;
ax.YAxis.TickLabels = {};
xlabel('Time (s)');
ylabel('Voltage (\muV)');

subplot(3,1,2);
plot(-0.800:0.001:0.199,beta(1001:2000),'LineWidth',thickPlotLineWidth,'Color',plotColours(2,:)); 
set_axis(gca,axisSettings);
ax = gca;
ax.YAxis.TickLabels = {};
xlabel('Time (s)');
ylabel('Voltage (\muV)');

subplot(3,1,3);
plot(beta(2001:3000),'LineWidth',thickPlotLineWidth,'Color',plotColours(3,:));
set_axis(gca,axisSettings);
ax = gca;
ax.YAxis.TickLabels = {};
% ax.XAxis.Visible = 'off';
% ax.YAxis.Visible = 'off';
% ax.YLim = commonYLim;
%title('Cue-to-Response (Scaled-Time)');
ax.XAxis.TickLabels = {'1','20','40','60','80','100'};
xlabel('Interval proportion (%)');
ylabel('Voltage (\muV)');

print(fullfile(tiffFolder,['fig_01c_recoveredcomps_' whichSimulation '.tiff']),'-dtiff',tiffRes);

%% Fig. XX (model third component as fixed rather than scaled)
gap = [0.16 0.12];
marg_h = [0.12 0.06];
marg_w = [0.1 0.06];
subplot = @(m,n,p) subtightplot (m, n, p,gap,marg_h, marg_w);
%   Input arguments (defaults exist):
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

thickPlotLineWidth = 1;
thinPlotLineWidth = 1; 
commonYLim = [-2 4];
commonWindow = 2800:4200;
makefigure(8,8);

subplot(3,1,1);
plot(0:0.001:0.999,betaFixedOnly(1:1000),'LineWidth',thickPlotLineWidth,'Color',plotColours(1,:)); 
set_axis(gca,axisSettings);
ax = gca;
ax.YAxis.TickLabels = {};
xlabel('Time (s)');
ylabel('Voltage (\muV)');

subplot(3,1,2);
plot(-0.800:0.001:0.199,betaFixedOnly(1001:2000),'LineWidth',thickPlotLineWidth,'Color',plotColours(2,:)); 
% xlim([-400 100]);
set_axis(gca,axisSettings);
ax = gca;
ax.YAxis.TickLabels = {};
xlabel('Time (s)');
ylabel('Voltage (\muV)');

subplot(3,1,3);
plot(0:0.001:2.999,betaFixedOnly(2001:5000),'LineWidth',thickPlotLineWidth,'Color',plotColours(3,:));
set_axis(gca,axisSettings);
ax = gca;
ax.YAxis.TickLabels = {};
xlabel('Time (s)');
ylabel('Voltage (\muV)');

print(fullfile(tiffFolder,['fig_01_recoveredcompsfixedonly_' whichSimulation '.tiff']),'-dtiff',tiffRes);

disp('model fits');
disp(['with scaled component: ', num2str(errorTS)]);
disp(['with only fixed components: ', num2str(errorFixedOnly)])

%% Fig. 1d
meanEpochedData = squeeze(mean(epoched_dat,2));

gap = [0.16 0.12];
marg_h = [0.12 0.04];
marg_w = [0.14 0.06];
subplot = @(m,n,p) subtightplot (m, n, p,gap,marg_h, marg_w);
%   Input arguments (defaults exist):
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 


fig = makefigure(8,6.5);

subplot(2,1,1);
for i = 1:length(intervals)
    plot((0:1000)./1000,squeeze(mean(epoched_dat(i,:,500:1500),2)),'LineWidth',plotLineWidth,'Color',eventColours(i,:));
    hold on;
end

xlabel('Time (s)');
ylabel('Voltage (\muV)');
l = legend('Short','Medium','Long');
l.Location = 'SouthEast';
l.Box = 'off';
set_axis(gca,axisSettings);
ax = gca;
ax.YLim = [-2 2];

subplot(2,1,2);
timesInt = [-800,200];
times = (timesInt(1):timesInt(2))./1000;
for i = 1:length(intervals)
    thisPlotTimes = intervals(i)-300 : intervals(i)+700;
    plot((timesInt(1):timesInt(2))./1000, squeeze(mean(epoched_dat(i,:,thisPlotTimes),2)),'LineWidth',plotLineWidth,'Color',eventColours(i,:));
    hold on;
end

% title('Response');
xlabel('Time (s)');
xlim(timesInt./1000);
ylabel('Voltage (\muV)');
l = legend('Short','Medium','Long');
l.Location = 'NorthWest';
l.Box = 'off';
set_axis(gca,axisSettings);

print(fullfile(tiffFolder,['fig_01d_simulatederps_' whichSimulation '.tiff']),'-dtiff',tiffRes);
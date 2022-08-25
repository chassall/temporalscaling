% Figure 4, Scale Factors
%
% Other m-files required: subtightplot, sigstar, notboxplot

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% 25-Aug-2022

close all; clear all; clc;

% Set results and output folders (change as needed)
resultsFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\results';
outputFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\figures\tiffs';

subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.25 0.15], [0.15 0.1], [0.1 0.01]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

tasks = {'production','perception','prediction','dm'};

allPs = [];
allTs = [];
allDFs = [];
allDs = [];
makefigure(9,7);
axs = {};
alsoAxs = [];
whichMeas = 'r2'; % 'r' or 'r2'
% Load scaling factor results for each task and plot
for taskI = 1:length(tasks)
    load(fullfile(resultsFolder,['scalefactor_' tasks{taskI} '.mat']),'allRPre','allRFixed','allRScaled','allR2Pre','allR2Fixed','allR2Scaled');
    
    switch whichMeas
        case 'r'
            sfPre = allRPre;
            sfFixed = allRFixed;
            sfScaled = allRScaled;
        case 'r2'
            sfPre = allR2Pre;
            sfFixed = allR2Fixed;
            sfScaled = allR2Scaled;
    end
    
    collPre = zeros(size(allRPre,1),1);
    collFixed = zeros(size(allRPre,1),1);
    collScaled = zeros(size(allRPre,1),1);
    switch tasks{taskI}
        case {'production','perception'} % short, medium, long conditions
            whichCombos = [3 1; 3 2];
        case {'prediction','dm'} % short and long conditions only
            whichCombos = [2 1];
    end

    for i = 1:size(whichCombos,1)
        collPre = collPre + sfPre(:,whichCombos(i,1),whichCombos(i,2));
        collFixed = collFixed + sfFixed(:,whichCombos(i,1),whichCombos(i,2));
        collScaled = collScaled + sfScaled(:,whichCombos(i,1),whichCombos(i,2));
    end
    collPre = collPre ./ size(whichCombos,1);
    collFixed = collFixed ./ size(whichCombos,1);
    collScaled = collScaled ./ size(whichCombos,1);

    
    subplot(2,2,taskI);
    thisNBP = notBoxPlot([collPre collFixed collScaled],'interval','tInterval');
    formatNBP(thisNBP);
    axs{taskI} = gca;
    alsoAxs(taskI) = gca;
    collDiff = collScaled - collFixed;
    [h,p,CI,stats] = ttest(collDiff);
    allPs(taskI) = p;
    allDFs(taskI)= stats.df;
    allTs(taskI) = stats.tstat;
    allDs(taskI) = mean(collDiff)/std(collDiff);
end

linkaxes(alsoAxs);

%% Formatting

plot_settings;
titles = {'Production','Perception','Prediction','Decision Making'};
for i = 1:length(axs)
    
    set(gcf,'CurrentAxes',axs{i});
    if allPs(i) < 0.05
    H = sigstar([2 3],allPs(i));
    set(H,'LineWidth',0.5);
    end
    % yline(0,':');
    
    axs{i}.Title.String = titles{i};
    axs{i}.Title.Position(2) = axs{i}.Title.Position(2) * 1.05;
    % axs{i}.Title.FontSize = plotFontSize; 
    axs{i}.Title.FontWeight = 'normal';
    axs{i}.Title.FontSize = plotFontSize;
    axs{i}.Box = 'off';
    axs{i}.FontSize = plotFontSize;
    axs{i}.XTickLabel = {'Original','Fixed-only','Scaled-only'};
    switch whichMeas
        case 'r'
            axs{i}.YLabel.String = 'Scaling Index (r)';
        case 'r2'
            axs{i}.YLabel.String = 'Scaling Index (R^2)';
    end
    axs{i}.XTickLabelRotation = 30;
    % axs{i}.YLim = [-0.8 1.4];
end

disp(allTs);
disp(allDFs);
disp(allPs);
disp(allDs);

print(fullfile(outputFolder,'fig_4.tiff'),'-dtiff','-r600');
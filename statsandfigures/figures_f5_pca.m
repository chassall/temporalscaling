% Figure 5, PCA and Residuals
%
% Other m-files required: subtightplot, notboxplot, plot_settings

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% 25-Aug-2022

%% Figures
close all; clear all;
% Set results and tiffs folders (change as needed)
resultsFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\results';
tiffFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\figures\tiffs';
tiffRes = '-r600';
load(fullfile(resultsFolder,'residual_results_production_task.mat'));
close all; % There may be some figures that were loaded

plot_settings;
short = squeeze(mean(doShortERP,1));
medium = squeeze(mean(doMediumERP,1));
long = squeeze(mean(doLongERP,1));
shortFixed = squeeze(mean(doShortERPF,1));
mediumFixed = squeeze(mean(doMediumERPF,1));
longFixed = squeeze(mean(doLongERPF,1));
shortScaled= squeeze(mean(doShortERPS,1));
mediumScaled = squeeze(mean(doMediumERPS,1));
longScaled = squeeze(mean(doLongERPS,1));

iElectrode = 15;
Electrode = 24;
figW = 5;
figH = 2.5;
axs = {};
legendStrings = {'Short','Medium','Long'};

shortTimes = doShortBeepInterval(1):(1/EEG.srate):(doShortBeepInterval(2)-(1/EEG.srate));
mediumTimes = doMediumBeepInterval(1):(1/EEG.srate):(doMediumBeepInterval(2)-(1/EEG.srate));
longTimes = doLongBeepInterval(1):(1/EEG.srate):(doLongBeepInterval(2)-(1/EEG.srate));

fig1 = makefigure(figW+2.5,figH);
plot(shortTimes,short(iElectrode,:),'LineWidth',plotLineWidth','Color',eventColours(1,:)); hold on;
plot(mediumTimes,medium(iElectrode,:),'LineWidth',plotLineWidth,'Color',eventColours(2,:));
plot(longTimes,long(iElectrode,:),'LineWidth',plotLineWidth,'Color',eventColours(5,:));
axs{1} = gca;
legend(legendStrings,'Box','off','Location','EastOutside');

fig2 = makefigure(figW,figH);
plot(shortTimes,shortFixed(iElectrode,:),'LineWidth',plotLineWidth,'Color',eventColours(1,:)); hold on;
plot(mediumTimes,mediumFixed(iElectrode,:),'LineWidth',plotLineWidth,'Color',eventColours(2,:));
plot(longTimes,longFixed(iElectrode,:),'LineWidth',plotLineWidth,'Color',eventColours(5,:));
axs{2} = gca;

fig3 = makefigure(figW,figH);
plot(shortTimes,shortScaled(iElectrode,:),'LineWidth',plotLineWidth,'Color',eventColours(1,:)); hold on;
plot(mediumTimes,mediumScaled(iElectrode,:),'LineWidth',plotLineWidth,'Color',eventColours(2,:));
plot(longTimes,longScaled(iElectrode,:),'LineWidth',plotLineWidth,'Color',eventColours(5,:));
axs{3} = gca;

linkaxes([axs{:}]);

for i = 1:length(axs)
    ax = axs{i};
    ax.FontSize = plotFontSize;
    ax.Box = 'off';
    ax.YLabel.String = 'Voltage (\muV)';
    ax.XLabel.String = 'Time (s)';
end

print(fig1, fullfile(tiffFolder,'fig_3_production_erp.tiff'),'-dtiff',tiffRes);
print(fig2, fullfile(tiffFolder,'fig_3_production_erp_fixed.tiff'),'-dtiff',tiffRes);
print(fig3, fullfile(tiffFolder,'fig_3_production_erp_scaled.tiff'),'-dtiff',tiffRes);

%% Load PCA results (production task)
close all;
clearvars -except resultsFolder tiffFolder tiffRes
load(fullfile(resultsFolder,'pca_results_production_task.mat'));

%% Mean Variance
tVal = abs(tinv(0.025,nParticipants-1));
meanVar = squeeze(mean(explainedByParticipant(:,:,1:2),1));
stdVar = squeeze(std(explainedByParticipant(:,:,1:2),[],1));
ciVar = tVal * stdVar ./ sqrt(nParticipants);
disp('Production Task: Mean Variance Explained (Supplementary Table 5');
disp(meanVar);
disp(meanVar - ciVar);
disp(meanVar + ciVar);

%% Plot PCA scores for the production task

figW = 5;
figH = 2.5;
plot_colours = cbrewer('seq','YlGnBu',nQuantiles+1);
makefigure(figW,figH);
hm = notBoxPlot(meanScores,'interval','tinterval');
formatNBP(hm,plot_colours);
ylabel('PC2 score (a.u.)');
xlabel('Response time quantile');
thisTable = array2table(meanScores,'VariableNames',quantileNames);
within = array2table(withinNames,'VariableNames',{'Quantile'});
rm = fitrm(thisTable,rmModel, 'WithinDesign',within);
[ranovatblb] = ranova(rm, 'WithinModel','Quantile');
disp('PCA Results: Production Task (Results Section)');
disp(ranovatblb{'(Intercept):Quantile','F'});
disp(ranovatblb{'(Intercept):Quantile','pValue'});

% Effect size
etaSquared = ranovatblb.SumSq(3) / sum(ranovatblb.SumSq);
partialEtaSquared = ranovatblb.SumSq(3) / (ranovatblb.SumSq(3) + ranovatblb.SumSq(4));
generalEtaSquared = ranovatblb.SumSq(3) / (ranovatblb.SumSq(3) + ranovatblb.SumSq(4) + ranovatblb.SumSq(2));
disp(partialEtaSquared);
disp(generalEtaSquared);

[H, pValue, SWstatistic] = swtest(meanScores(:,1));
[H, pValue, SWstatistic] = swtest(meanScores(:,2));
[H, pValue, SWstatistic] = swtest(meanScores(:,3));

ax = gca;
set_axis(ax,axisSettings);
% ax.XTickLabel = {'very fast','fast','on time','slow','very slow'};
ax.XTickLabel = {'early','"on time"','late'};
print(fullfile(tiffFolder,['fig_5_' whichTask '_' whichAnalysis '_pca_pc2_by_condition.tiff']),'-dtiff',tiffRes);

%% Mean PC2 Scores
disp('Mean PC2 Scores: Production Task (Supplementary Table 6');
meanData = mean(meanScores);
tVal = abs(tinv(0.025,nParticipants-1));
stdData = std(meanScores);
ciData = tVal * stdData ./ sqrt(nParticipants);
disp(meanData);
disp(meanData - ciData);
disp(meanData + ciData);

%% Find and plot mean PCA scores for each quantile (for reviewer 1)
meanPCAScores = squeeze(mean(allPCAScores(:,:,whichChannel,:),1));
stdPCAScores = squeeze(std(allPCAScores(:,:,whichChannel,:),[],1));
tVal = abs(tinv(0.025,nParticipants-1));
ciPCAScpres = tVal * stdPCAScores ./ sqrt(nParticipants);

participantPCAScores = squeeze(mean(allPCAScores(:,:,whichChannel,:),1));

plot_colours = cbrewer('seq','YlGnBu',nQuantiles+1);
h = {};

subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.14], [0.18 0.14], [0.08 0.02]);
fig = makefigure(19,7);

for c = 1:3
    disp(conditionStrings{c});
    subplot(1,3,c);
    theseScores = squeeze(allPCAScores(:,c,whichChannel,:));
    h{c} = notBoxPlot(theseScores);
    formatNBP(h{c},plot_colours);
    title(conditionStrings{c});
    ylabel('PC2 score (a.u.)');
    xlabel('Response time quantile');
    
end

print(fullfile(tiffFolder,['fig_5_' whichTask '_' whichAnalysis '_pca_by_condition_and_quantile.tiff']),'-dtiff',tiffRes);

%% Big plot of effect of +/- PC2 (Supplement)
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.10 0.08], [0.1 0.06], [0.08 0.02]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis.
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins
%

poolParticipants = 0;

makefigure(17.17,9);

plot_settings;
whichSlots = [1 2 3];

for c = 1:3
    subplot(2,3,whichSlots(c))
    
    if poolParticipants
        % Plot coeff from pooled data
        thisPlot = plot(trialTimes{c},allCoeff{c}(:,1:2),'LineWidth',plotLineWidth);
        for pcaI = 1:2
            thisPlot(pcaI).Color = plotColours(pcaI,:);
        end
        
    else
        % Plot mean coeff (unpooled)
        smoothPCACoeff{c} = movmean(allPCACoeff{c},mmeanwindow,2);
        meanPCACoeff = squeeze(mean(smoothPCACoeff{c},1));
        stdPCACoeff = squeeze(std(smoothPCACoeff{c},[],1));
        ciPCACoeff = tVal*stdPCACoeff./sqrt(nParticipants);
        b1 = boundedline(trialTimes{c},meanPCACoeff(:,1),ciPCACoeff(:,1),'cmap',plotColours(1,:), 'alpha');
        b1.LineWidth = plotLineWidth;
        hold on;
        b2 = boundedline(trialTimes{c},meanPCACoeff(:,2),ciPCACoeff(:,2),'cmap',plotColours(2,:), 'alpha');
        b2.LineWidth = plotLineWidth;
    end
    
    ax = gca;
    set_axis(ax, axisSettings);
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'Comp. amp. (a.u.)';
    ax.XLim = [0 2.5];
    if poolParticipants
        l1 = legend('PC1','PC2');
    else
        l1 = legend([b1,b2],'PC1','PC2');
    end
    l1.Box = 'off';
    l1.Location = 'best';
    
    if ismember(c,[1 2 3])
        ax.XLabel.String = '';
    end
    
    if ismember(c,[2 3])
        ax.YLabel.String = '';
    end
end

meanScores = squeeze(mean(allPCAScores(:,:,whichChannel,:),2));

% Effect of PCA2
whichPlots = [4 5 6];
plot_settings;
pc2ScoresToPlot = [];
allMinMax = [];
for c = 1:3
    if poolParticipants
        thisMin = -25;
        thisMax = 25;
    else
        thisMin = round(mean(min(allPCAScores2{c}(:,:,2),[],2)));
        thisMax = round(mean(max(allPCAScores2{c}(:,:,2),[],2)));
        allMinMax(c,1) = thisMin;
        allMinMax(c,2) = thisMax;
    end
    pc2ScoresToPlot = [pc2ScoresToPlot; linspace(thisMin,thisMax,7)];
end
disp('Min/Max PC2 scores for the figure');
disp(allMinMax);

for c = 1:3
    subplot(2,3,whichPlots(c));
    plotColours = cbrewer('div','RdYlBu',size(pc2ScoresToPlot,2));
    plotColours = cbrewer('div','RdBu',size(pc2ScoresToPlot,2));
    plotColours = flip(plotColours);
    
    plotHandles = [];
    meanPC1Shape = mean(theseScores(:,1)*allCoeff{c}(:,1)',1);
    
    % Indiviual PCAs
    meanPCACoeffInd = squeeze(mean(smoothPCACoeff{c},1));
    stdPCACoeffInd = squeeze(std(smoothPCACoeff{c},[],1));
    ciPCACoeffInd = tVal*stdPCACoeffInd./sqrt(nParticipants);
    
    meanPC1ShapeInd = [];
    meanPC2CoeffInd = [];
    for iParticipant = 1:nParticipants
        thisPC1Coeff = smoothPCACoeff{c}(iParticipant,:,1);
        thisPC1Score = squeeze(allPCAScores2{c}(iParticipant,:,1));
        thisPC1Rec = thisPC1Score' * thisPC1Coeff;
        thisPC1Rec = mean(thisPC1Rec,1);
        meanPC1ShapeInd = [meanPC1ShapeInd; thisPC1Rec];
        
        thisPC2Coeff = smoothPCACoeff{c}(iParticipant,:,2);
        meanPC2CoeffInd = [meanPC2CoeffInd; thisPC2Coeff];
    end
    grandMeanPC1ShapeInd = mean(meanPC1ShapeInd,1);
    grandMeanPC2CoeffInd = mean(meanPC2CoeffInd,1);
    
    for w = 1:size(pc2ScoresToPlot,2)
        if poolParticipants
            thisPC2Shape = mean(pc2ScoresToPlot(c,w)*allCoeff{c}(:,2)',1);
            thisData = meanPC1Shape + thisPC2Shape;
        else
            % Individual
            thisPC2Shape = pc2ScoresToPlot(c,w)*grandMeanPC2CoeffInd;
            thisData = grandMeanPC1ShapeInd + thisPC2Shape;
            
        end
        
        [~,thisMin] = min(thisData);
        plotHandles(w) = plot(trialTimes{c},thisData,'Color',plotColours(w,:),'LineWidth',plotLineWidth);
        hold on;
        xlabel('Time (s)');
        ylabel('Voltage (\muV)');
        ax = gca;
        set_axis(ax, axisSettings);
        if poolParticipants
            % ax.YLim(1) = ax.YLim(1) - 0.5;
            % ax.YLim = [-1 0.5];
        else
            ax.YLim  = [-15 2];
        end
        ax.XLim = [0 2.5];
    end
    if poolParticipants
        plotHandles(size(pc2ScoresToPlot,2)+1) = plot(trialTimes{c},meanPC1Shape,'k','LineWidth',plotLineWidth);
    else
        % Individual
        plotHandles(size(pc2ScoresToPlot,2)+1) = plot(trialTimes{c},grandMeanPC1ShapeInd,'k','LineWidth',plotLineWidth);
    end
    
    if ismember(c,[2 3])
        ax.YLabel.String = '';
    end
end
print(fullfile(tiffFolder,['sfig_4' whichTask '_' whichAnalysis '_add_subtract_pc2.tiff']),'-dtiff',tiffRes);

%% Load PCA scores for the DM task
close all;
clearvars -except resultsFolder tiffFolder tiffRes

load(fullfile(resultsFolder,'dm_pca.mat'),'quantileData','nParticipants','explainedByParticipant');

figW = 5;
figH = 2.5;
makefigure(figW,figH);
plot_colours = cbrewer('seq','YlGnBu',3);
plot_settings;
hm = notBoxPlot(quantileData,'interval','tinterval');
formatNBP(hm, plot_colours);
ylabel('PC2 score (a.u.)');
xlabel('Decision time quantile');
ax = gca;
set_axis(ax,axisSettings);
% ax.XTickLabel = {'very fast','fast','on time','slow','very slow'};
ax.XTickLabel = {'early','average','late'};
print(fullfile(tiffFolder,'fig_5_dm_pca.tiff'),'-dtiff',tiffRes);

disp('Mean PC2 Scores: Decision-Making (Supplementary Table 6');
meanData = mean(quantileData);
tVal = abs(tinv(0.025,nParticipants-1)); 
stdData = std(quantileData);
ciData = tVal * stdData ./ sqrt(nParticipants);
disp(meanData);
disp(meanData - ciData);
disp(meanData + ciData);

%% Stats
thisTable = array2table(quantileData,'VariableNames',{'q1','q2','q3'});
within = array2table({'1';'2';'3'},'VariableNames',{'Quantile'});
rm = fitrm(thisTable,'q1-q3 ~ 1', 'WithinDesign',within);
[ranovatblb] = ranova(rm, 'WithinModel','Quantile');
disp('PCA Results: DM Task (Results Section)');
disp(ranovatblb{'(Intercept):Quantile','F'});
disp(ranovatblb{'(Intercept):Quantile','pValue'});

% Effect size
etaSquared = ranovatblb.SumSq(3) / sum(ranovatblb.SumSq);
partialEtaSquared = ranovatblb.SumSq(3) / (ranovatblb.SumSq(3) + ranovatblb.SumSq(4));
generalEtaSquared = ranovatblb.SumSq(3) / (ranovatblb.SumSq(3) + ranovatblb.SumSq(4) + ranovatblb.SumSq(2));
disp(partialEtaSquared);
disp(generalEtaSquared);

[H, pValue, SWstatistic] = swtest(quantileData(:,1));
[H, pValue, SWstatistic] = swtest(quantileData(:,2));
[H, pValue, SWstatistic] = swtest(quantileData(:,3));

%% Mean Variance
tVal = abs(tinv(0.025,nParticipants-1));
meanVar = squeeze(mean(explainedByParticipant(:,1:2),1));
stdVar = squeeze(std(explainedByParticipant(:,1:2),[],1));
ciVar = tVal * stdVar ./ sqrt(nParticipants);
disp('Mean Variance Explained: Decision-Making Task (Supplementary Figure 5)');
disp(meanVar);
disp(meanVar - ciVar);
disp(meanVar + ciVar);
% Supplementary Figure 1, Different Interpolation Methods
%
% Other m-files required: subtightplot

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% 25-Aug-2022

close all; clear all; clc;

% Set output folder, change as needed
tiffFolder = 'C:\Users\chassall\OneDrive - Nexus365\Projects\2020_EEG_TimeEstimationTS_Hassall\figures\tiffs';

widthRERP = 20;
testX = eye(widthRERP);
methods = {'box','bilinear','nearest'};
allTogether = {};
plot_settings;
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.15 0.1], [0.1 0.1]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

myColourMap = cbrewer('seq','Greys',256);
for iMethod = 1:length(methods)
    
    method = methods{iMethod};
    reallyShortX = imresize(testX,[widthRERP*0.25,widthRERP],'Method',method);
    shorterX = imresize(testX,[widthRERP*0.5,widthRERP],'Method',method);
    longerX = imresize(testX,[widthRERP*1.5,widthRERP],'Method',method);
    reallyLongX = imresize(testX,[widthRERP*1.75,widthRERP],'Method',method);
    
    f = make_figure(13.6,5.5);
    % figure('Name',method);
    ax(1) = subplot(3,6,[1,2]);
    imagesc(shorterX,[0,1]); title('Compressed'); xlabel('\beta index'); ylabel('EEG time (ms)');
    ax(2) = subplot(3,6,[3,4,9,10]);
    imagesc(testX,[0,1]); title('Original Stick Function'); xlabel('\beta index'); ylabel('EEG time (ms)');
    ax(3) = subplot(3,6,[5,6,11,12,17,18]);
    imagesc(longerX,[0,1]); title('Stretched'); xlabel('\beta index'); ylabel('EEG time (ms)');
    colormap(myColourMap);
    
     ax(4) = subplot(3,6,[13,14]);
     cb = colorbar();
     cb.Location = 'south';
     colormap(myColourMap);
     ax(4).XAxis.Visible = 'off';
     ax(4).YAxis.Visible = 'off';
     cb.Ticks = [0,0.5,1];
     cb.Label.String = 'Regressor value';

    for i = 1:3
        set_axis(ax(i),axisSettings);
    end
    % print(method,'-dtiff','-r300');
    
    allTogether{iMethod} = [reallyShortX; shorterX; testX; longerX; reallyLongX];
end

numColours = 255;
midColour = ceil(numColours/2);
% myColourMap = cbrewer('div','RdBu',numColours);
myColourMap = cbrewer('seq','Greys',256);
% myColourMap = flip(myColourMap);
% myColourMap(midColour,:) = [1 1 1];
clear ax;
make_figure(14.94,11.81);
make_figure(10.94,11.81);
for i = 1:length(methods)
    
    ax(i) = subplot(1,length(methods),i);
    imagesc(allTogether{i},[0,1]); title(methods{i}); colormap(myColourMap);
    xlabel('\beta index'); ylabel('EEG time (ms)');
end
for i = 1:length(ax)
    set_axis(ax(i),axisSettings);
end
print(fullfile(tiffFolder,'fig_sf1_stretchmethods.tiff'),'-dtiff','-r600');

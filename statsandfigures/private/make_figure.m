function fig = make_figure(width,height)
%MAKE_FIGURE Make a publication-ready figure
% Width (in cm) can be 3.0, 9.0, 14.0, or 19.0 aka minimum size, single column, 1.5-column, or double column
% Height can be 24.0 cm max
% Note that wide PowerPoint slides are 33.87 wide by 19.05 high
% Note that Nature defines single-column as 8.9 cm, double as 18.3 cm, and
% 1.5 column as 12.0-13.6 cm
% 

fig = figure;
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.Units = 'centimeters';
padding = 0.1;
fig.PaperPosition = [0, 0, width, height];
fig.PaperSize = [width, height];
fig.Position = [padding, padding, width-padding, height-padding];
fig.InvertHardcopy = 'off';
fig.Color = 'white';

end


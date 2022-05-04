function [prettyFig, text_color] = format3DPlots(fig, mode, plottitle, ...
    xAxisLabel, yAxisLabel, zAxisLabel)
%FORMATPLOTS Stephen Heirtzler's line plot formatter.
%   Applies nice formatting to a line plot.

% Pick between dark and light mode
if ~exist("mode","var") || mode ~= "dark" && mode ~= "light"
    mode = "dark";
end

set(fig, 'InvertHardcopy', 'off')
ax = fig.CurrentAxes;

% Switch statement to apply mode colors
switch mode
    case "dark"
        text_color = "#f9f9f9"; ax_color = "#f2f2f2";
        fig.Color = "#00080b"; ax.Color = "#00080b";
        ax.XColor = ax_color; ax.YColor = ax_color; ax.ZColor = ax_color;
    case "light"
        text_color = "#00080b"; ax_color = "#080808";
        fig.Color = "#ffffff"; ax.Color = "#f5f5f5";
        ax.XColor = ax_color; ax.YColor = ax_color; ax.ZColor = ax_color;
end

ax.TickLabelInterpreter = 'latex';
ax.FontSize = 9;
ax.XDir = 'reverse';
% Add title and axis labels if present
if exist("plottitle","var")
    title(ax, plottitle,'Interpreter','latex','Color',text_color,'FontSize',9)
end
if exist("xAxisLabel","var")
    xlabel(ax, xAxisLabel,'Interpreter','latex','Color',text_color)
end
if exist("yAxisLabel","var")
    ylabel(ax, yAxisLabel,'Interpreter','latex','Color',text_color)
end
if exist("zAxisLabel","var")
    zlabel(ax, zAxisLabel,'Interpreter','latex','Color',text_color)
end


% Remove unnecessary clutter from figure
fig.MenuBar = 'none';
fig.ToolBar = 'none';
grid on
view(45,30)
% Get Screen dimensions to center figure in screen
set(0,'units','inches')
Inch_SS = get(0,'screensize');

% Resize figure to fit in a standard sized two column document
fig.Units = 'inches';
fig.Position = [(Inch_SS(3)/2 - (3.13/2)) (Inch_SS(4)/2 - (2.17/2)) 3.13 2.17];

prettyFig = fig;
end

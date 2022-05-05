function plotstress(normal_x, normal_y, shear_xy)
%PLOTSTRESS plots stress components in the style of the textbook
fig = figure;
ax = gca;

set(fig, 'InvertHardcopy', 'off')
axis equal
ax.XLim = [0 1]; ax.YLim = [0 1];
C = permute(ones(4,3).*linspace(0.8,0.7,4)',[1 3 2]);
patch([0.25 0.75 0.75 0.25],[0.75 0.75 0.25 0.25],C,'LineWidth',1.2)
ax.Position = [0 0 1 1];
axis off
fig.Color = "#ffffff"; ax.Color = "#ffffff";

prefixes = ["","k","M","G"];
powers =  [1,1e3,1e6,1e9];

% Arrows
%Shear
if abs(shear_xy) > 1
annotation("arrow", ...
    [0.5+0.1.*sign(shear_xy+eps) 0.5+0.3.*sign(shear_xy+eps)],[0.82 0.82], ...
    'HeadLength',20,'HeadStyle','plain','LineWidth',2,'Color','#e1373d')
annotation("arrow", ...
    [0.5-0.1.*sign(shear_xy+eps) 0.5-0.3.*sign(shear_xy+eps)],[0.18 0.18], ...
    'HeadLength',20,'HeadStyle','plain','LineWidth',2,'Color','#e1373d')

annotation("arrow", ...
    [0.18 0.18],[0.5-0.1.*sign(shear_xy+eps) 0.5-0.3.*sign(shear_xy+eps)], ...
    'HeadLength',20,'HeadStyle','plain','LineWidth',2,'Color','#e1373d')
annotation("arrow", ...
    [0.82 0.82],[0.5+0.1.*sign(shear_xy+eps) 0.5+0.3.*sign(shear_xy+eps)], ...
    'HeadLength',20,'HeadStyle','plain','LineWidth',2,'Color','#e1373d')
end

%Normal
if abs(normal_x) > 1
annotation("arrow", ...
    [0.85-0.1.*sign(normal_x+eps) 0.85+0.1.*sign(normal_x+eps)],[0.5 0.5], ...
    'HeadLength',20,'HeadStyle','plain','LineWidth',2,'Color','#e1373d')
annotation("arrow", ...
    [0.15+0.1.*sign(normal_x+eps) 0.15-0.1.*sign(normal_x+eps)],[0.5 0.5], ...
    'HeadLength',20,'HeadStyle','plain','LineWidth',2,'Color','#e1373d')
end

if abs(normal_y) > 1
annotation("arrow", ...
    [0.5 0.5],[0.85-0.1.*sign(normal_y+eps) 0.85+0.1.*sign(normal_y+eps)], ...
    'HeadLength',20,'HeadStyle','plain','LineWidth',2,'Color','#e1373d')
annotation("arrow", ...
    [0.5 0.5],[0.15+0.1.*sign(normal_y+eps) 0.15-0.1.*sign(normal_y+eps)], ...
    'HeadLength',20,'HeadStyle','plain','LineWidth',2,'Color','#e1373d')
end

% Labels
%Shear
if abs(shear_xy) > 1
sheartext = sprintf('%0.1f %sPa', ...
    shear_xy/powers(floor(log10(abs(shear_xy))/3 + 1)), ...
    prefixes(floor(log10(abs(shear_xy))/3 + 1)));
text(0.5+0.2.*sign(shear_xy+eps),0.85,sheartext, ...
    'Interpreter','latex', 'FontSize',12,'Color','#e1373d', ...
    'HorizontalAlignment','center','VerticalAlignment','bottom')
end

%Normal
if abs(normal_y) > 1
normytext = sprintf('%0.1f %sPa', ...
    normal_y/powers(floor(log10(abs(normal_y))/3 + 1)), ...
    prefixes(floor(log10(abs(normal_y))/3 + 1)));
text(0.5-0.12.*sign(shear_xy+eps),0.85,normytext, ...
    'Interpreter','latex', 'FontSize',12,'Color','#e1373d', ...
    'HorizontalAlignment','center','VerticalAlignment','bottom')
end

if abs(normal_x) > 1
normxtext = sprintf('%0.1f %sPa', ...
    normal_x/powers(floor(log10(abs(normal_x))/3 + 1)), ...
    prefixes(floor(log10(abs(normal_x))/3 + 1)));
text(0.87,0.5-0.07.*sign(shear_xy+eps),normxtext, ...
    'Interpreter','latex', 'FontSize',12,'Color','#e1373d', ...
    'HorizontalAlignment','center','VerticalAlignment','middle')
end

set(0,'units','inches')
Inch_SS = get(0,'screensize');

% Remove unnecessary clutter from figure
%fig.MenuBar = 'none';
fig.ToolBar = 'none';

% Resize figure to fit in a standard sized two column document
fig.Units = 'inches';
fig.Position = [(Inch_SS(3)/2 - (3/2)) (Inch_SS(4)/2 - (3/2)) 3 3];

annotation("arrow",[0.5 0.5],[0.5 0.65],'HeadLength',10,'LineWidth',1)
annotation("arrow",[0.5 0.65],[0.5 0.5],'HeadLength',10,'LineWidth',1)
text(0.52,0.65,'$y\prime$','Interpreter','latex')
text(0.65,0.53,'$x\prime$','Interpreter','latex')

end


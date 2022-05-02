%% Calculate Positions (Run this before plotting positions)
clear;clc
% WARNING: REQUIRES SYMBOLIC MATH TOOLBOX

% All lengths are in meters
% All angles are in radians

group_num = 17;
% Any variable with s appended is a solution array for every crankangle
phi_2s = linspace(0,2*pi,360);
w_2 = deg2rad(group_num*5 + 20);
alpha_2 = 0;

% Intialize known lengths and angles
r_5 = 0.2;
r_2 = 0.3;
r_4 = 1.4;
r_bo = 0.75;
r_y = 0.5;
phi_bo = deg2rad(360-106);

% Use symbolic math toolbox to initialize unknown position variables
syms r_4ac phi_4 r_x phi_5 phi_2 real
assume(phi_5 < 0);
assumeAlso(phi_5 > -sym(pi/2))

% Vector loop A_o - A - C - D, expressed in polar form
equ1 = r_2*exp(1i*phi_2) + r_4ac*exp(1i*phi_4) == r_x + r_y*1i;
% The second equation is the complex conjugate of the first, both equations
% are valid, and therefore can be used to solve the system.
equ1con=r_2*exp(-1i*phi_2) + r_4ac*exp(-1i*phi_4) == r_x + r_y*-1i;

% Vector loop A_o - B_o - B - C - D, also expressed in polar form
equ2 = r_bo*exp(1i*phi_bo) + r_5*exp(1i*phi_5) + r_4*exp(1i*phi_4) == r_x + r_y*1i;
% Same story here
equ2con=r_bo*exp(-1i*phi_bo) + r_5*exp(-1i*phi_5) + r_4*exp(-1i*phi_4) == r_x + r_y*-1i;

% Initialize empty arrays for the unknown angles.
r_4acs = zeros(size(phi_2s));
phi_4s = zeros(size(phi_2s));
phi_5s = zeros(size(phi_2s));

% Iterate through every value of phi_2 and evaluate the system of equations
% at each angle. This finds the positions.
W = waitbar(0,'Calulating Positions'); % User feedback :)
for n = 1:numel(phi_2s)
    % Substitute the value of phi_2, then numerically solve the system with
    % vpasolve.
    P = vpasolve(subs(equ1,phi_2,phi_2s(n)), ...
                 subs(equ1con,phi_2,phi_2s(n)) ...
                 ,equ2,equ2con, ...
        r_4ac, phi_4, r_x, phi_5, [r_4/2; pi/2; 0; -pi/8]);
    % Write unknown angles and lengths to arrays
    r_4acs(n) = P.r_4ac;
    phi_4s(n) = P.phi_4;
    phi_5s(n) = P.phi_5;
    if mod(n,3)
        waitbar(n/numel(phi_2s),W);
    end
end
clear equ1 equ1con equ2 equ2con
delete(W)

r_2x = r_2.*sin(phi_2s);
r_2y = r_2.*cos(phi_2s);

r_5x = r_5.*sin(phi_5s);
r_5y = r_5.*cos(phi_5s);

r_4x = r_4.*sin(phi_4s);
r_4y = r_4.*cos(phi_4s);
%% Calulate Velocities (run this before plotting velocities)
clc
% MAKE SURE TO RUN THE SECTION ABOVE FIRST
try phi_4s(1); catch %#ok<NOEFF> 
    error('No position data. Please run Calulate Positions section first.')
end
% Use symbolic math toolbox to initialize unknown velocity variables
syms w_4 w_5 v_x v_4ac real

equ3 = r_2*w_2*exp(1i*(phi_2 + pi/2)) + r_4ac*w_4*exp(1i*(phi_4 + pi/2)) ...
         + v_4ac*exp(1i*phi_4) == v_x;
equ3con = r_2*w_2*exp(-1i*(phi_2 + pi/2)) + r_4ac*w_4*exp(-1i*(phi_4 + pi/2)) ...
         + v_4ac*exp(-1i*phi_4) == v_x;

equ4 = r_5*w_5*exp(1i*(phi_5 + pi/2)) + r_4*w_4*exp(1i*(phi_4 + pi/2)) == v_x;
equ4con = r_5*w_5*exp(-1i*(phi_5 + pi/2)) + r_4*w_4*exp(-1i*(phi_4 + pi/2)) == v_x;

w_4s = zeros(size(phi_2s));
w_5s = zeros(size(phi_2s));
v_xs = zeros(size(phi_2s));
v_4acs = zeros(size(phi_2s));

% Iterate through every value of phi_2. Then evaluate the
% system of equations at each angle. This finds the velocities.
W = waitbar(0,'Calulating Velocities'); % User feedback :)
for n = 1:numel(phi_2s)
    % Substitute the value of phi_2, then numerically solve the system with
    % vpasolve.
    V = vpasolve(subs(equ3,[phi_2 phi_4 r_4ac],[phi_2s(n) phi_4s(n) r_4acs(n)]), ...
                 subs(equ3con,[phi_2 phi_4 r_4ac],[phi_2s(n) phi_4s(n) r_4acs(n)]), ...
                 subs(equ4,[phi_4 phi_5],[phi_4s(n) phi_5s(n)]), ...
                 subs(equ4con,[phi_4 phi_5],[phi_4s(n) phi_5s(n)]), ...
                 w_4, w_5, v_x, v_4ac);
    w_4s(n) = V.w_4;
    w_5s(n) = V.w_5;
    v_xs(n) = V.v_x;
    v_4acs(n) = V.v_4ac;
    if mod(n,3)
        waitbar(n/numel(phi_2s),W);
    end
end
clear equ3 equ3con equ4 equ4con
delete(W)
%% Calulate Accelerations (run this before plotting accelerations)
clc
% MAKE SURE TO RUN THE SECTION ABOVE FIRST
try w_4s(1); catch %#ok<NOEFF> 
    error('No velocity data. Please run Calulate Velocities section first.')
end
% Use symbolic math toolbox to initialize unknown acceleration variables
syms alpha_4 alpha_5 a_4ac a_x real

equ5 = r_2*w_2^2*exp(1i*(phi_2+pi)) + (2*v_4ac*w_4 + r_4ac*alpha_4)*exp(1i*(phi_4+pi/2)) ...
    + a_4ac*exp(1i*phi_4) + r_4ac*w_4^2*exp(1i*(phi_4+pi)) == a_x;

equ5con = r_2*w_2^2*exp(-1i*(phi_2+pi)) + (2*v_4ac*w_4 + r_4ac*alpha_4)*exp(-1i*(phi_4+pi/2)) ...
    + a_4ac*exp(-1i*phi_4) + r_4ac*w_4^2*exp(-1i*(phi_4+pi)) == a_x;

equ6 = r_5*alpha_5*exp(1i*(phi_5+pi/2)) + r_5*w_5^2*exp(1i*(phi_5+pi)) ...
     + r_4*alpha_4*exp(1i*(phi_4+pi/2)) + r_4*w_4^2*exp(1i*(phi_4+pi)) == a_x;

equ6con = r_5*alpha_5*exp(-1i*(phi_5+pi/2)) + r_5*w_5^2*exp(-1i*(phi_5+pi)) ...
     + r_4*alpha_4*exp(-1i*(phi_4+pi/2)) + r_4*w_4^2*exp(-1i*(phi_4+pi)) == a_x;

alpha_4s = zeros(size(phi_2s));
alpha_5s = zeros(size(phi_2s));
a_xs = zeros(size(phi_2s));
a_4acs = zeros(size(phi_2s));

% Iterate through every value of phi_2. Then evaluate the
% system of equations at each angle. This finds the velocities.
W = waitbar(0,'Calulating Accelerations'); % User feedback :)
for n = 1:numel(phi_2s)
    % Substitute the value of phi_2, then numerically solve the system with
    % vpasolve.
    A = vpasolve(subs(equ5,[phi_2 phi_4 r_4ac v_4ac w_4], ...
                           [phi_2s(n) phi_4s(n) r_4acs(n) v_4acs(n) w_4s(n)]), ...
                 subs(equ5con,[phi_2 phi_4 r_4ac v_4ac w_4], ...
                              [phi_2s(n) phi_4s(n) r_4acs(n) v_4acs(n) w_4s(n)]), ...
                 subs(equ6,[phi_4 phi_5 w_4 w_5], ...
                           [phi_4s(n) phi_5s(n) w_4s(n) w_5s(n)]), ...
                 subs(equ6con,[phi_4 phi_5 w_4 w_5], ...
                              [phi_4s(n) phi_5s(n) w_4s(n) w_5s(n)]), ...
                 alpha_4, alpha_5, a_x, a_4ac);
    
    alpha_4s(n) = A.alpha_4;
    alpha_5s(n) = A.alpha_5;
    a_xs(n) = A.a_x;
    a_4acs(n) = A.a_4ac;
    if mod(n,3)
        waitbar(n/numel(phi_2s),W);
    end
end
clear equ6 equ6con equ5 equ5con
delete(W)

l_4=linspace(0,r_4,49);     %[m]
l_2=linspace(0,r_2,49);     %[m]
l_5=linspace(0,r_5,49);     %[m]
%%
[L_2,PHI_2s] = meshgrid(l_2,phi_2s);
[L_4,PHI_4s] = meshgrid(l_4,phi_4s);
[L_5,PHI_5s] = meshgrid(l_5,phi_5s);
[~,ALPHA_5s] = meshgrid(l_2,real(alpha_5s));
[~,ALPHA_4s] = meshgrid(l_2,real(alpha_4s));
[~,W_4s] = meshgrid(l_2,real(w_4s));
[~,W_5s] = meshgrid(l_2,real(w_5s));

A_2s = (L_2/2).* w_2.^2 .* exp(1i.*(PHI_2s + pi));
A_5s = (L_5/2).*W_5s.^2 .* exp(1i.*(PHI_5s + pi)) ...
    + ALPHA_5s.*(L_5/2).* exp(1i.*(PHI_5s+pi/2));
% a_4acs = a_5.*2 + ((r_4ac/2).*w_4s.^2+a_4acs) .* exp(1i.*(phi_4s+pi)) ...
%     + alpha_4s.*(r_4ac/2).* exp(1i.*(phi_4s+pi/2));
A_4s = A_5s.*2 + (L_4/2).*W_4s.^2 .* exp(1i.*(PHI_4s+pi)) ...
    + ALPHA_4s.*(L_4/2).* exp(1i.*(PHI_4s+pi/2));

A_2x = real(A_2s);
A_2y = imag(A_2s);

A_5x = real(A_5s);
A_5y = imag(A_5s);

A_4x = real(A_4s);
A_4y = real(A_4s);




%% Dynamic Force Calculations


p = 7800; %kg/m^3
d = .1; %m
a = 0.05; %m
b = .2; %m

m_2 = p*(pi*((d/2)^2)*r_2);%kg
m_4 = p*(a*b*r_4); %kg
m_5 = p*(pi*(d/2)^2*r_5); %kg
m_6 = 20; %kg

m5 = 1; %N/m
R = 1000; %N



f_6x = m_6*A_4x + R;
f_6y = 0;
f_6 = vecnorm(cat(3,f_6x, zeros(size(f_6x))),2,3);

f_4x = R - m_5.*A_5x + m_4.*A_4x;
f_4y = m_2.*A_2y - m_4.*A_4y;
f_4 = vecnorm(cat(3,f_4x, f_4y),2,3);

f_2x = f_4x - f_4x + m_2.*A_2x;
f_2y = f_4y - f_4y + m_2.*A_2y;
f_2 = vecnorm(cat(3,f_2x, f_2y),2,3);

f_3x = f_4x-f_2x;
f_3y = f_2y-f_4y;
f_3 = vecnorm(cat(3,f_3x, f_3y),2,3);

f_5x = f_4y - f_4y + m_5.*A_5x;
f_5y = f_4y - f_4y + m_5.*A_5y;
f_5 = vecnorm(cat(3,f_5x, f_5y),2,3);

[~,R_2x] = meshgrid(l_2,r_2x);
[~,R_2y] = meshgrid(l_2,r_2y);
[~,R_4x] = meshgrid(l_2,r_4x);
[~,R_4y] = meshgrid(l_2,r_4y);
[~,R_5x] = meshgrid(l_2,r_5x);
[~,R_5y] = meshgrid(l_2,r_5y);

M2 = -f_2x.*(R_2y./2)+f_2y.*(R_2x./2)+f_2x.*(R_2y./2)-f_2y.*(R_2x./2);
M3 = m5+f_5x.*(R_4y./2)+f_5y.*(R_4x./2)+R.*(R_4y./2);
M4 = m5+f_5x.*(R_4y./2)+f_5y.*(R_4x./2)+R.*(R_4y./2);
M5 = m5;
M6 = 0;

%% Position Plots
close all; clc
% Reassemble Solved angles and lengths into matrix of positions vs. time
membArray = [r_2.*exp(1i.*phi_2s); cumsum([r_bo*exp(1i*phi_bo).*ones(1,numel(phi_2s)); ...
    r_5.*exp(1i.*phi_5s); ...
    r_4.*exp(1i.*phi_4s); ...
    ],1)];

% Draw Trajectory lines
fig = figure;
trA = line(real(membArray(1,1:end)),imag(membArray(1,1:end)));
trB = line(real(membArray(3,1:end)),imag(membArray(3,1:end)));
trC = line(real(membArray(4,1:end)),imag(membArray(4,1:end)));
trA.LineWidth = 1.5; trA.Color = [0.0 0.3 0.5].*1;
trA.DisplayName = "A";
trB.LineWidth = 1.5; trB.Color = [0.6 0.3 0.0].*1.4;
trB.DisplayName = "B";
trC.LineWidth = 1.5; trC.Color = [0.5 0.4 0.0].*1.6;
trC.DisplayName = "C";

formatPlots(fig,'light',[],"$x$ [m]","$y$ [m]");

axis equal
axis padded
xticks(-1:0.5:1);
legend("Location","southeast","Interpreter","latex")
fig.Position(1:2) = fig.Position(1:2) - fig.Position(3:4).*0.6;
fig.Position(3:4) = fig.Position(3:4).*2;

% Draw position component for point C
fig2 = figure;
ax = gca;
trC = line(phi_2s.*180/pi,real(membArray(4,1:end)));
trC.LineWidth = 1.5;
trC.Color = [0.8 0.4 0.5].*0.7;

formatPlots(fig2,'light',[],"$\theta$ [degrees]","$x$ [m]");

ax.XLimitMethod = "tight";
ax.YLimitMethod = "tickaligned";
xticks(0:45:360);
fig2.Position(1:2) = fig2.Position(1:2) - fig2.Position(3:4).*0.5;
fig2.Position(3:4) = fig2.Position(3:4).*2;
%% Velocity Plot
close all; clc
fig3 = figure;
ax = gca;
velC = line(phi_2s.*180/pi,real(v_xs));
velC.LineWidth = 1.5;
velC.Color = [0.1 0.3 0.6].*1.2;
formatPlots(fig3,'light',[],"$\theta$ [degrees]","$v_x \; [\mathrm{m}/\mathrm{s}]$");
ax.XLimitMethod = "tight";
ax.YLimitMethod = "tickaligned";
xticks(0:45:360);
fig3.Position(1:2) = fig3.Position(1:2) - fig3.Position(3:4).*0.5;
fig3.Position(3:4) = fig3.Position(3:4).*2;
%% Acceleration Plot
close all; clc

a_2s = (r_2/2).* w_2.^2 .* exp(1i.*(phi_2s+pi));
a_5s = (r_5/2).*w_5s.^2 .* exp(1i.*(phi_5s+pi)) ...
    + alpha_5s.*(r_5/2).* exp(1i.*(phi_5s+pi/2));

a_4s = a_5s.*2 + (r_4/2).*w_4s.^2 .* exp(1i.*(phi_4s+pi)) ...
    + alpha_4s.*(r_4/2).* exp(1i.*(phi_4s+pi/2));

fig4 = figure;
ax = gca;
acc2 = line(phi_2s.*180/pi,abs(a_2s));
acc4 = line(phi_2s.*180/pi,abs(a_4s));
acc5 = line(phi_2s.*180/pi,abs(a_5s));
acc6 = line(phi_2s.*180/pi,abs(a_xs));

acc2.LineWidth = 1.5; acc2.Color = [0.0 0.3 0.5];
acc2.DisplayName = "2";
acc4.LineWidth = 1.5; acc4.Color = [0.1 0.7 0.6];
acc4.DisplayName = "4";
acc5.LineWidth = 1.5; acc5.Color = [0.9 0.45 0.0];
acc5.DisplayName = "5";
acc6.LineWidth = 1.5; acc6.Color = [0.5 0.4 0.0].*1.5;
acc6.DisplayName = "6";
formatPlots(gcf,'light',[],"$\theta$ [degrees]","$$|\vec{a}| \; [\mathrm{m}/\mathrm{s}^2]$$");
ax.XLimitMethod = "tight";
ax.YLimitMethod = "tickaligned";
xticks(0:45:360);
legend("Location","northwest","Interpreter","latex")
fig4.Position(1:2) = fig4.Position(1:2) - fig4.Position(3:4).*0.5;
fig4.Position(3:4) = fig4.Position(3:4).*2;

%% Forces Plotting

figure()
plot(phi_2s.*180/pi,f_2(:,25)',"LineWidth",1.5)
hold on
plot(phi_2s.*180/pi,f_3(:,25)',"LineWidth",1.5)
plot(phi_2s.*180/pi,f_4(:,25)',"LineWidth",1.5)
plot(phi_2s.*180/pi,f_5(:,25)',"LineWidth",1.5)
plot(phi_2s.*180/pi,f_6(:,25)',"LineWidth",1.5)
hold off
formatPlots(gcf,'light')
xlim([0 360])
xticks(0:45:360);
ylim([-200 1300])
legend('Force on Member 2','Force on Member 3','Force on Member 4','Force on Member 5','Force on Member 6','location','east')

print('-dpng','-r300')
%% Axial stress plotting
ax = cos(PHI_2s);
ay = sin(PHI_2s);
a = sqrt(ax.^2 + ay.^2);

ax1 = cos(PHI_4s);
ay1 = sin(PHI_4s);
b = sqrt(ax1.^2 + ay1.^2);

ax2 = cos(PHI_5s);
ay2 = sin(PHI_5s);
c = sqrt(ax2.^2 + ay2.^2);

F_2a = f_2.*a;
F_4a = f_4.*b;
F_5a = f_5.*c;

figure('Name','Axial_Force_2');
ax = gca;
surf(L_2,PHI_2s.*180/pi,real(F_2a),"LineWidth",1.5,"EdgeColor","none");
format3DPlots(gcf,'light','Axial force on member 2','Member length','Crank Angle','Axial force')
ylim([0 360])
yticks(0:90:360);
print('-dpng','-r300')

figure('Name','Axial_Force_4');
surf(L_4,PHI_2s.*180/pi,real(F_4a),"LineWidth",1.5,"EdgeColor","none");
format3DPlots(gcf,'light','Axial force on member 4','Member length','Crank Angle','Axial force')
ylim([0 360])
yticks(0:90:360);
print('-dpng','-r300')

figure('Name','Axial_Force_5');
surf(L_5,PHI_2s.*180/pi,real(F_5a),"LineWidth",1.5,"EdgeColor","none");
format3DPlots(gcf,'light','Axial force on member 5','Member length','Crank Angle','Axial force')
ylim([0 360])
yticks(0:90:360);
print('-dpng','-r300')

%% Shear plotting

axi = cos((PHI_2s+(pi/2)));
ayi = sin((PHI_2s+(pi/2)));
w = sqrt(axi.^2 + ayi.^2);


axj = cos(PHI_4s+(pi/2));
ayj = sin(PHI_4s+(pi/2));
x = sqrt(axj.^2 + ayj.^2);

axk = cos(PHI_5s+(pi/2));
ayk = sin(PHI_5s+(pi/2));
y = sqrt(axk.^2 + ay.^2);


F_2b = f_2.*w;
F_4b = f_4.*x;
F_5b = f_5.*y;

figure('Name','Shear_Force_2');
surf(L_2,PHI_2s.*180/pi,real(F_2b),"LineWidth",1.5,"EdgeColor","none")
format3DPlots(gcf,'light','Shear force on member 2','Member length','Crank Angle','Shear force')
ylim([0 360])
yticks(0:120:360);
%print('-dpng','-r300')

figure('Name','Shear_Force_4');
surf(L_4,PHI_2s.*180/pi,real(F_4b),"LineWidth",1.5,"EdgeColor","none")
format3DPlots(gcf,'light','Shear force on member 4','Member length','Crank Angle','Shear force')
ylim([0 360])
yticks(0:120:360);
%print('-dpng','-r300')

figure('Name','Shear_Force_5');
surf(L_5,PHI_2s.*180/pi,real(F_5b),"LineWidth",1.5,"EdgeColor","none")
format3DPlots(gcf,'light','Shear force on member 5','Member length','Crank Angle','Shear force')
ylim([0 360])
yticks(0:120:360);
%print('-dpng','-r300')
%% Bending Plotting
% Cameron Charette

m5=linspace(0,1,49);

[M5,~] = meshgrid(m5,phi_2s);
figure;
surf(L_2,PHI_2s.*180/pi,M2,"EdgeColor","none")
format3DPlots(gcf,'light','Internal Bending Moment for Member 2','Member length','Crank Angle');
zlim([0 0.5])
clim([0 0.5])
ylim([0 360])
yticks(0:120:360);
print('-dpng','-r300')
figure;
surf(L_4,PHI_2s.*180/pi,M4,"EdgeColor","none")
format3DPlots(gcf,'light','Internal Bending Moment for Member 4','Member length','Crank Angle');
ylim([0 360])
yticks(0:120:360);
print('-dpng','-r300')
figure;
surf(L_5,PHI_2s.*180/pi,M5,"EdgeColor","none")
format3DPlots(gcf,'light','Internal Bending Moment for Member 5','Member length','Crank Angle');
ylim([0 360])
yticks(0:120:360);
print('-dpng','-r300')

%% Functions
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


function [prettyFig, text_color] = formatPlots(fig, mode, plottitle, ...
    xAxisLabel, yAxisLabel)
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
        fig.Color = "#ffffff"; ax.Color = "#ffffff";
        ax.XColor = ax_color; ax.YColor = ax_color; ax.ZColor = ax_color;
end

ax.TickLabelInterpreter = 'latex';
ax.FontSize = 9;

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

% Remove unnecessary clutter from figure
fig.MenuBar = 'none';
fig.ToolBar = 'none';

% Get Screen dimensions to center figure in screen
set(0,'units','inches')
Inch_SS = get(0,'screensize');

% Resize figure to fit in a standard sized two column document
fig.Units = 'inches';
fig.Position = [(Inch_SS(3)/2 - (3.13/2)) (Inch_SS(4)/2 - (2.17/2)) 3.13 2.17];

prettyFig = fig;
end
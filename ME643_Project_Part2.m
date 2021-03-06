clear;clc;close all % Clear enviroment (clear all is slow and unnecessary)
% Grab all the data variables from deliverable 1.
% This saves an enormous amount of time.
load('Deliverable_1_data.mat')

%% FATIGUE ANALYSIS
% BEAM CROSS SECTIONS
a = 0.05;    % [m] height
b = 0.2;     % [m] width
d2 = d;      % [m] beam 2 diameter
d5 = d;      % [m] beam 5 diameter
Area_4 = a*b;
Area_2 = pi*(d2/2)^2;
Area_5 = pi*(d5/2)^2;

% STRESS CALCULATIONS
% Beam 2
F_2 = F_2a + M2;
F_2avg = mean(F_2,1);

idx2 = find(F_2avg == max(F_2avg));
S_2 = F_2(:,idx2) ./ Area_2; % Stress on beam 2 at the highest stress location

% Beam 4
F_4 = F_4a + M4;
F_4avg = mean(F_4,1);

idx4 = find(F_4avg == max(F_4avg));
S_4 = F_4(:,idx4) ./ Area_4; % Stress on beam 4 at the highest stress location

% Beam 5
F_5 = F_5a + M5;
F_5avg = mean(F_5,1);

idx5 = find(F_5avg == max(F_5avg));
S_5 = F_5(:,idx5) ./ Area_5; % Stress on beam 5 at the highest stress location

% PLOT Stress vs Phi for highest stress point of each member
fig1 = figure(1);
ax = gca;
line(phi_2s*(180/pi),S_2);
formatPlots(fig1,'light',[],"$\theta$ [degrees]","$\sigma$ [Pa]");
ylim([0 max(1200)])
ax.XLimitMethod = "tight";
xticks(0:45:360);

fig2 = figure(2);
ax = gca;
line(phi_2s*(180/pi),S_4)
formatPlots(fig2,'light',[],"$\theta$ [degrees]","$\sigma$ [Pa]");
ax.XLimitMethod = "tight"; ax.YLimitMethod = "tickaligned";
xticks(0:45:360);

fig3 = figure(3);
ax = gca;
line(phi_2s*(180/pi),S_5)
formatPlots(fig3,'light',[],"$\theta$ [degrees]","$\sigma$ [Pa]");
ax.XLimitMethod = "tight"; ax.YLimitMethod = "tickaligned";
xticks(0:45:360);

% CALCULATE alternating stress (S_a) and mean stress (S_m)

% Beam 2
S_2a = 0.5*(max(S_2) - min(S_2));
S_2m = 0.5*(max(S_2) + min(S_2));

% Beam 4
S_4a = 0.5*(max(S_4) - min(S_4));
S_4m = 0.5*(max(S_4) + min(S_4));

% Beam 5
S_5a = 0.5*(max(S_5) - min(S_5));
S_5m = 0.5*(max(S_5) + min(S_5));

%% Deflections
%Uses code for deliverable 1
E= 200000e3; %Modulus of elasticity 
L2 = 0.3; %[m]
L4 = 1.4; %[m]
L5 = 0.2; %[m]

r2 = d2/2; %radius of member 2
r5 = d5/2; %radius of member 5

I2 = (pi*r2^(4))/4; %moment of inertia member 2
I4 = (a*b^(3))/12; %moment of inertia member 4
I5 = (pi*r5^(4))/4; %moment of inertia member 5

%Deflection for member 2
y2max = (m_2*L2^(2))/ (9*sqrt(3)*E*I2);  %Deflection on member 2

% %P= center load on beam
P4 = m_2;
% %Deflection for member 4
y4max = (P4*L4^(3))/(48*E*I4);  %Deflection on memeber 4

%Deflection for member 5
y5max = (m_5*L5^(2))/ (9*sqrt(3)*E*I5);  %Deflection on member 5

ycr2 = L2/360;
ycr4 = L4/360; 
ycr5 = L5/360; 

N_deflect_2 = ycr2/y2max;
N_deflect_4 = ycr4/y4max;
N_deflect_5 = ycr5/y5max;

%% Buckling Analysis
% Stephen Heirtzler
E= 200000e3;

% Use eular column buckling formula on each column
I_2 = (pi*d2^4)/64;
P_cr_buckling_2 = (pi^2*E*I_2)/(r_2^2);
N_buckling_2 = P_cr_buckling_2/max(mean(F_2a,2),[],"all");

I_4 = (b*a^3)/12;
P_cr_buckling_4 = (pi^2*E*I_4)/(r_4^2);
N_buckling_4 = P_cr_buckling_4/max(mean(F_4a,2),[],"all");

I_5 = (pi*d5^4)/64;
P_cr_buckling_5 = (pi^2*E*I_5)/(r_5^2);
N_buckling_5 = P_cr_buckling_5/max(mean(F_5a,2),[],"all");

%% Pin tearout and shear analysis
% Stephen Heirtzler
%Bearing
d_pin = 0.02;

A_bearing = pi/4 * d * d_pin;

S_max_bearing = (max(F_4(:,end))/A_bearing);
N_bearing = 86e6 / S_max_bearing;

%Tearout
A_tear = ((b-d_pin)*a)^2;
S_max_tearout = (max(F_4(:,end))/A_tear);
N_tearout = 86e6 / S_max_tearout;

%Direct shear
A_shear = pi*(d_pin/2)^2;
S_max_shear = (max(F_4(:,end)*0.5)/A_shear);
N_shear = 86e6 / S_max_shear;

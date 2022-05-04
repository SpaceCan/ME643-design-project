close all
clear all;clc
% Grab all the data variables from deliverable 1.
% This saves an enormous amount of time.
load('Deliverable_1_data.mat')

%% FATIGUE ANALYSIS
% BEAM CROSS SECTIONS
a = 0.05; % [m] beam 4 width
b = 0.05; % [m] beam 4 thickness
d2 = 0.05; % [m] beam 2 diameter
d5 = 0.05; % [m] beam 2 diameter
Area_4 = a*b;
Area_2 = pi*(d2/2)^2;
Area_5 = pi*(d5/2)^2;

% STRESS CALCULATIONS
% Beam 2
for i = 1:length(f_2(:,1))
    for j = 1:length(f_2(1,:))
        F_2(i,j) = F_2a(i,j) + M2(i,j);
        F_2avg(j) = mean(F_2(:,j));
    end
end
idx2 = find(F_2avg == max(F_2avg));
S_2 = F_2(:,idx2) ./ Area_2; % Stress on beam 2 at the highest stress loacation
% Beam 4
for i = 1:length(f_4(:,1))
    for j = 1:length(f_4(1,:))
        F_4(i,j) = F_4a(i,j) + M4(i,j);
        F_4avg(j) = mean(F_4(:,j));
    end
end
idx4 = find(F_4avg == max(F_4avg));
S_4 = F_4(:,idx4) ./ Area_4; % Stress on beam 4 at the highest stress loacation
% Beam 5
for i = 1:length(f_5(:,1))
    for j = 1:length(f_5(1,:))
        F_5(i,j) = F_5a(i,j) + M5(i,j);
        F_5avg(j) = mean(F_5(:,j));
    end
end
idx5 = find(F_5avg == max(F_5avg));
S_5 = F_5(:,idx5) ./ Area_5; % Stress on beam 5 at the highest stress loacation

% PLOT Stress vs Phi for highest stress point of each member
figure
plot(phi_2s*(180/pi),S_2)
figure
plot(phi_2s*(180/pi),S_4)
figure
plot(phi_2s*(180/pi),S_5)

% CALCULATE alternating stress (S_a) and mean stress (S_m)
% Beam 2
S_2max = max(S_2);
S_2min = min(S_2);
S_2a = 0.5*(S_2max - S_2min);
S_2m = 0.5*(S_2max + S_2min);
% Beam 4
S_4max = max(S_4);
S_4min = min(S_4);
S_4a = 0.5*(S_4max - S_4min);
S_4m = 0.5*(S_4max + S_4min);
% Beam 5
S_5max = max(S_5);
S_5min = min(S_5);
S_5a = 0.5*(S_5max - S_5min);
S_5m = 0.5*(S_5max + S_5min);



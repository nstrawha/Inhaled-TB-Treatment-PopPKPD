%% To recreate Figure 5 in main text drug concentrations in lung
clear all
% Parameters the same for all 4 drugs (from Tables S2–S5, S7)
BW = 70; % kg
phys = loadPhysiology(BW);
kF = 0.252; % gut transit rate 

%Rifampicin 
ka_RIF = 1.08; % absorption rate [1/h] - supplement says 
CL_RIF = 7.86; % systemic clearance [L/h] - supplement says 
fR_RIF = 0.1830; % fractional renal clearance - supplement says 
kr_RIF = 0.17; % gut reabsorption rate [1/h]
pt_RIF = loadPartitionCoefficients('rifampicin');

% Ethambutol
ka_EMB = 0.22; % absorption rate [1/h]
CL_EMB = 49.93; % systemic clearance [L/h] - supplement says 49.99
fR_EMB = 0.79; % fractional renal clearance
kr_EMB = 0; % gut reabsorption rate [1/h]
pt_EMB = loadPartitionCoefficients('ethambutol');

% Isoniazid
% Slow
ka_INHS = 4.11; % absorption rate [1/h] 4.11 slow; 2.86 fast
CL_INHS = 9.17; % systemic clearance [L/h] supplement says - 9.16 slow; 
fR_INHS = 0.29; % fractional renal clearance 0.07 fast; 0.29 slow
kr_INHS = 0; % gut reabsorption rate [1/h] 
% Fast 
ka_INHF = 2.89; % absorption rate [1/h] 4.11 slow; supplement says - 2.86 fast
CL_INHF = 24.34; % systemic clearance [L/h] 9.16 slow; supplement says - 24.56 fast
fR_INHF = 0.07; % fractional renal clearance 0.07 fast; 0.29 slow
kr_INHF = 0; % gut reabsorption rate [1/h] 
pt_INH = loadPartitionCoefficients('isoniazid');

% Pyrazinamide
ka_PYZ = 1.39; % absorption rate [1/h] - supplement says 1.36
CL_PYZ = 4.14; % systemic clearance [L/h] - supplement says 4.1
fR_PYZ = 0.09; % fractional renal clearance
kr_PYZ = 0; % gut reabsorption rate [1/h]
pt_PYZ = loadPartitionCoefficients('pyrazinamide');


options = odeset('RelTol',1e-6,'AbsTol',1e-8); %To solve ODEs

%% Figure 5  (lung)  Day 1 for RIF, INH, PYZ and Day 7 for EMB
fig = figure();
fig.Position = [00 00 1920 1080];
tiledlayout(2,2);
set(0,'DefaultFigureWindowStyle','docked');


% RIF - Day 1 
A0_RIF = zeros(1,18); A0_RIF(18) = 600;
[t_RIF, A_RIF] = ode23s(@(t, A) humanEqns(t, A, ka_RIF, kr_RIF,kF, CL_RIF, fR_RIF, phys, pt_RIF), [0:0.01:30], A0_RIF,options);
nexttile
plot(t_RIF, A_RIF(:,3),'LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (\mug/mL)');
title('Rifampicin');
set(gca,'FontSize',25);
grid on;

% EMB - Day 7 
nexttile
A0_EMB= zeros(1,18); A0_EMB(18) = 1200;
for d = 1:7
    [t_EMB, A_EMB] = ode23s(@(t, A) humanEqns(t, A, ka_EMB, kr_EMB,kF, CL_EMB, fR_EMB, phys, pt_EMB), [0 24], A0_EMB,options);
    A0_EMB=[A_EMB(end,1:17)';A_EMB(end,18)+1200];
end
plot(t_EMB, A_EMB(:,3),'LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (\mug/mL)');
title('Ethambutol');
set(gca,'FontSize',25);
grid on;

%INH - Day 1 step is 0.1 instead of 0.01 
A0_INH = zeros(1,18); A0_INH(18) = 300;
[t_INHF, A_INHF] = ode23s(@(t, A) humanEqns(t, A, ka_INHF, kr_INHF,kF, CL_INHF, fR_INHF, phys, pt_INH), [0:0.1:24], A0_INH,options);
[t_INHS, A_INHS] = ode23s(@(t, A) humanEqns(t, A, ka_INHS, kr_INHS,kF, CL_INHS, fR_INHS, phys, pt_INH), [0:0.1:24], A0_INH,options);

nexttile
plot(t_INHF, A_INHF(:,3), 'DisplayName', 'Fast','LineWidth',2); hold on
plot(t_INHS, A_INHS(:,3), 'DisplayName', 'Slow','LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (\mug/mL)');
title('Isoniazid');
set(gca,'FontSize',25);
legend show;
grid on;

% PYZ - Day 1 
A0_PYZ= zeros(1,18); A0_PYZ(18) = 1500;
[t_PYZ, A_PYZ] = ode23s(@(t, A) humanEqns(t, A, ka_PYZ, kr_PYZ,kF, CL_PYZ, fR_PYZ, phys, pt_PYZ), [0 48], A0_PYZ,options);
nexttile
plot(t_PYZ, A_PYZ(:,3),'LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (\mug/mL)');
ylim( [0 100])
title('Pyrazinamide'); 
set(gca,'FontSize',25);
grid on;
sgtitle('Lung concentration','FontSize',30);
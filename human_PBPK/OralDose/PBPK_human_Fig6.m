%% To recreate Figure 6 in main text drug concentrations at various extrapulmonary sites
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

%% Figure 6 (2 parts, 4 panels then 1 solo panel)
fig = figure();
fig.Position = [00 00 1920 1080];
tiledlayout(2,2);
set(0,'DefaultFigureWindowStyle','docked');

% RIF
A0_RIF = zeros(1,18); A0_RIF(18) = 600;
for d = 1:7
    [t_RIF, A_RIF] = ode23s(@(t, A) humanEqns(t, A, ka_RIF, kr_RIF,kF, CL_RIF, fR_RIF, phys, pt_RIF), [0 24], A0_RIF,options);
    A0_RIF=[A_RIF(end,1:17)';A_RIF(end,18)+600];
   

end

nexttile
plot(t_RIF, A_RIF(:,16), 'DisplayName', 'Lymph Node','LineWidth',2); hold on;
plot(t_RIF, A_RIF(:,5), 'DisplayName', 'Brain','LineWidth',2);
plot(t_RIF, A_RIF(:,11), 'DisplayName', 'Bone','LineWidth',2);
plot(t_RIF, A_RIF(:,9), 'DisplayName', 'Skin','LineWidth',2);
plot(t_RIF, A_RIF(:,13), 'DisplayName', 'Kidney','LineWidth',2);
plot(t_RIF, A_RIF(:,15), 'DisplayName', 'Liver','LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (\mug/mL)');
title('Rifampicin');
set(gca,'FontSize',25);
legend show;
grid on;

% EMB
A0_EMB= zeros(1,18); A0_EMB(18) = 1200;
for d = 1:7
    [t_EMB, A_EMB] = ode23s(@(t, A) humanEqns(t, A, ka_EMB, kr_EMB,kF, CL_EMB, fR_EMB, phys, pt_EMB), [0 24], A0_EMB,options);
    A0_EMB=[A_EMB(end,1:17)';A_EMB(end,18)+1200];
end

nexttile
plot(t_EMB, A_EMB(:,16), 'DisplayName', 'Lymph Node','LineWidth',2); hold on;
plot(t_EMB, A_EMB(:,5), 'DisplayName', 'Brain','LineWidth',2);
plot(t_EMB, A_EMB(:,11), 'DisplayName', 'Bone','LineWidth',2);
plot(t_EMB, A_EMB(:,9), 'DisplayName', 'Skin','LineWidth',2);
plot(t_EMB, A_EMB(:,13), 'DisplayName', 'Kidney','LineWidth',2);
plot(t_EMB, A_EMB(:,15), 'DisplayName', 'Liver','LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (\mug/mL)');
title('Ethambutol');
set(gca,'FontSize',25);
legend show;
grid on;

%INH
A0_INHF = zeros(1,18); A0_INH(18) = 300;
for d = 1:7
    [t_INHF, A_INHF] = ode23s(@(t, A) humanEqns(t, A, ka_INHF, kr_INHF,kF, CL_INHF, fR_INHF, phys, pt_INH), [0 24], A0_INHF,options);
    A0_INHF=[A_INHF(end,1:17)';A_INHF(end,18)+300];
end

nexttile
plot(t_INHF, A_INHF(:,16), 'DisplayName', 'Lymph Node','LineWidth',2); hold on;
plot(t_INHF, A_INHF(:,5), 'DisplayName', 'Brain','LineWidth',2);
plot(t_INHF, A_INHF(:,11), 'DisplayName', 'Bone','LineWidth',2);
plot(t_INHF, A_INHF(:,9), 'DisplayName', 'Skin','LineWidth',2);
plot(t_INHF, A_INHF(:,13), 'DisplayName', 'Kidney','LineWidth',2);
plot(t_INHF, A_INHF(:,15), 'DisplayName', 'Liver','LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (\mug/mL)');
title('Isoniazid Fast');
set(gca,'FontSize',25);
legend show;
grid on;
A0_INHS = zeros(1,18); A0_INH(18) = 300;
for d = 1:7
    [t_INHS, A_INHS] = ode23s(@(t, A) humanEqns(t, A, ka_INHS, kr_INHS,kF, CL_INHS, fR_INHS, phys, pt_INH), [0 24], A0_INHS,options);
    A0_INHS=[A_INHS(end,1:17)';A_INHS(end,18)+300];
end

nexttile
plot(t_INHS, A_INHS(:,16), 'DisplayName', 'Lymph Node','LineWidth',2); hold on;
plot(t_INHS, A_INHS(:,5), 'DisplayName', 'Brain','LineWidth',2);
plot(t_INHS, A_INHS(:,11), 'DisplayName', 'Bone','LineWidth',2);
plot(t_INHS, A_INHS(:,9), 'DisplayName', 'Skin','LineWidth',2);
plot(t_INHS, A_INHS(:,13), 'DisplayName', 'Kidney','LineWidth',2);
plot(t_INHS, A_INHS(:,15), 'DisplayName', 'Liver','LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (\mug/mL)');
title('Isoniazid Slow');
set(gca,'FontSize',25);
legend show;
grid on;
sgtitle('Selected Tissue Conc Day 7','FontSize',30)

% PYZ
A0_PYZ= zeros(1,18); A0_PYZ(18) = 1600;
for d = 1:7
    [t_PYZ, A_PYZ] = ode23s(@(t, A) humanEqns(t, A, ka_PYZ, kr_PYZ,kF, CL_PYZ, fR_PYZ, phys, pt_PYZ), [0 24], A0_PYZ,options);
    A0_PYZ=[A_PYZ(end,1:17)';A_PYZ(end,18)+1600];
end


fig = figure();
fig.Position = [00 00 1920 1080];
set(0,'DefaultFigureWindowStyle','docked');

plot(t_PYZ, A_PYZ(:,16), 'DisplayName', 'Lymph Node','LineWidth',2); hold on;
plot(t_PYZ, A_PYZ(:,5), 'DisplayName', 'Brain','LineWidth',2);
plot(t_PYZ, A_PYZ(:,11), 'DisplayName', 'Bone','LineWidth',2);
plot(t_PYZ, A_PYZ(:,9), 'DisplayName', 'Skin','LineWidth',2);
plot(t_PYZ, A_PYZ(:,13), 'DisplayName', 'Kidney','LineWidth',2);
plot(t_PYZ, A_PYZ(:,15), 'DisplayName', 'Liver','LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (\mug/mL)');
ylim( [0 100])
title('Pyrazinamide - selected tissues');
set(gca,'FontSize',25);
legend show;
grid on;



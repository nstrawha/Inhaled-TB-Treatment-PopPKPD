
%% Option 1 - use partition coefficients from a mix of sources
BW = 0.025; % bodyweight of mouse in kg
phys = loadPhysiologyMouse(BW);
KP_1 = loadPartitionCoefficients('rifampicin');

%Rifampicin
fup = 0.029;    % Fraction unbound in plasma (Supp Table 1)
CLr = 19;       % Fraction of renal clearance fraction (Supp Table 1)
F = 1;          % Drug bioavailability (Supp Table 1)
ka = 0.3713;    % Rate of absorption (1/h) (Supp Table 3)It says the 'best fitted parameter'
CL = 0.037;     % Total body clearance (L/h) (Supp Table 3) It says the 'best fitted parameter'
BP = 0.9;       % Blood:plasma ratio

A0 = zeros(1,9); A0(7) = 10*BW; % drug is 10 mg/kg
[t_1, A_1] = ode15s(@(t, A) mouseEqns(t, A, fup, CLr,F, ka, CL,BP, phys,KP_1), [0 24], A0);

% Plot various drug concentrations
figure;

subplot(2,2,1);

% lung compartments
plot(t_1, A_1(:,1), 'DisplayName', 'Arterial','LineWidth',2); hold on;
plot(t_1, A_1(:,2), 'DisplayName', 'Venous','LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (mg/L)');
title('Plasma Concentration');
set(gca,'FontSize',25);
legend show;
grid on;

% arterial and venous concentrations
subplot(2,2,2);
plot(t_1, A_1(:,3), 'DisplayName', 'Lung Tissue','LineWidth',2); hold on;
xlabel('Time (h)'); ylabel('Concentration (mg/L)');
title('Drug Concentration in Lung Tissue');
set(gca,'FontSize',25);
legend show;
grid on;

subplot(2,2,3);
% oral and gut compartments
plot(t_1, A_1(:,7), 'DisplayName', 'Absorbed Drug (AD)','LineWidth',2);hold on;
plot(t_1, A_1(:,6), 'DisplayName', 'Gut','LineWidth',2);
xlabel('Time (h)'); ylabel('Amount of Drug (mg)');
title('Drug Amount in Dose Compartments');
set(gca,'FontSize',25);
legend show;
grid on;

sgtitle('Rifampicin - Mouse Day 1 - Model 1 Mix ','FontSize',30);

%% Option 2 - use the partition coefficients from Lyon 2013 for mouse

BW = 0.025; % bodyweight of mouse in kg
phys = loadPhysiologyMouse(BW);

%Rifampicin
fup = 0.029;    % Fraction unbound in plasma (Supp Table 1)
CLr = 19;       % Fraction of renal clearance fraction (Supp Table 1)
F = 1;          % Drug bioavailability (Supp Table 1)
ka = 0.3713;    % Rate of absorption (1/h) (Supp Table 3)It says the 'best fitted parameter'
CL = 0.037;     % Total body clearance (L/h) (Supp Table 3) It says the 'best fitted parameter'
BP = 0.9;       % Blood:plasma ratio


A0 = zeros(1,9); A0(7) = 10*BW; % drug is 10 mg/kg
KP_2 = loadPartitionCoefficients_Lyons('rifampin');
[t_2, A_2] = ode15s(@(t, A) mouseEqns(t, A, fup, CLr,F, ka, CL,BP, phys,KP_2), [0 24], A0);

% Plot various drug concentrations
figure;

subplot(2,2,1);

% lung compartments
plot(t_2, A_2(:,1), 'DisplayName', 'Arterial','LineWidth',2); hold on;
plot(t_2, A_2(:,2), 'DisplayName', 'Venous','LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (mg/L)');
title('Plasma Concentration');
set(gca,'FontSize',25);
legend show;
grid on;

% arterial and venous concentrations
subplot(2,2,2);
plot(t_2, A_2(:,3), 'DisplayName', 'Lung Tissue','LineWidth',2); hold on;
xlabel('Time (h)'); ylabel('Concentration (mg/L)');
title('Drug Concentration in Lung Tissue');
set(gca,'FontSize',25);
legend show;
grid on;

subplot(2,2,3);
% oral and gut compartments
plot(t_2, A_2(:,7), 'DisplayName', 'Absorbed Drug (AD)','LineWidth',2);hold on;
plot(t_2, A_2(:,6), 'DisplayName', 'Gut','LineWidth',2);
xlabel('Time (h)'); ylabel('Amount of Drug (mg)');
title('Drug Amount in Dose Compartments');
set(gca,'FontSize',25);
legend show;
grid on;

sgtitle('Rifampicin - Mouse Day 1 - Model 2 Lyons','FontSize',30);

%% Option 3 - try to calculate the partition coefficients using RR equations with info from Ramachandran
B = Blood_Properties(); %Blood cell properties
TC = RatTissueComposition(); %Rat Tissue Composition

%%%Drug Properties%%%
RIF = RIF_Properties(B); %Properties of Rifampicin

%%%Calculation of Partition coefficients%%%
KP_3 = calculateTissuePartition(RIF,TC,B); 

BW = 0.025; % bodyweight of mouse in kg
phys = loadPhysiologyMouse(BW);

%Rifampicin
fup = 0.029;    % Fraction unbound in plasma (Supp Table 1)
CLr = 19;       % Fraction of renal clearance fraction (Supp Table 1)
F = 1;          % Drug bioavailability (Supp Table 1)
ka = 0.3713;    % Rate of absorption (1/h) (Supp Table 3)It says the 'best fitted parameter'
CL = 0.037;     % Total body clearance (L/h) (Supp Table 3) It says the 'best fitted parameter'
BP = 0.9;       % Blood:plasma ratio


A0 = zeros(1,9); A0(7) = 10*BW; % drug is 10 mg/kg
[t_3, A_3] = ode15s(@(t, A) mouseEqns(t, A, fup, CLr,F, ka, CL,BP, phys,KP_3), [0 24], A0);

% Plot various drug concentrations
figure;

subplot(2,2,1);

% lung compartments
plot(t_3, A_3(:,1), 'DisplayName', 'Arterial','LineWidth',2); hold on;
plot(t_3, A_3(:,2), 'DisplayName', 'Venous','LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (mg/L)');
title('Plasma Concentration');
set(gca,'FontSize',25);
legend show;
grid on;

% arterial and venous concentrations
subplot(2,2,2);
plot(t_3, A_3(:,3), 'DisplayName', 'Lung Tissue','LineWidth',2); hold on;
xlabel('Time (h)'); ylabel('Concentration (mg/L)');
title('Drug Concentration in Lung Tissue');
set(gca,'FontSize',25);
legend show;
grid on;

subplot(2,2,3);
% oral and gut compartments
plot(t_3, A_3(:,7), 'DisplayName', 'Absorbed Drug (AD)','LineWidth',2);hold on;
plot(t_3, A_3(:,6), 'DisplayName', 'Gut','LineWidth',2);
xlabel('Time (h)'); ylabel('Amount of Drug (mg)');
title('Drug Amount in Dose Compartments');
set(gca,'FontSize',25);
legend show;
grid on;

sgtitle('Rifampicin - Mouse Day 1 - Model 3 Calculate KP','FontSize',30);
%% Option 4 - directly use partition coefficients from Ramachandran

BW = 0.025; % bodyweight of mouse in kg
phys = loadPhysiologyMouse(BW);

%Rifampicin
fup = 0.029;    % Fraction unbound in plasma (Supp Table 1)
CLr = 19;       % Fraction of renal clearance fraction (Supp Table 1)
F = 1;          % Drug bioavailability (Supp Table 1)
ka = 0.3713;    % Rate of absorption (1/h) (Supp Table 3)It says the 'best fitted parameter'
CL = 0.037;     % Total body clearance (L/h) (Supp Table 3) It says the 'best fitted parameter'
BP = 0.9;       % Blood:plasma ratio

A0 = zeros(1,9); A0(7) = 10*BW; % drug is 10 mg/kg
KP_4 =loadPartitionCoefficients_Ramachandran('rifampicin');
[t_4, A_4] = ode15s(@(t, A) mouseEqns(t, A, fup, CLr,F, ka, CL,BP, phys,KP_4), [0 24], A0);

% Plot various drug concentrations
figure;

subplot(2,2,1);

% lung compartments
plot(t_4, A_4(:,1), 'DisplayName', 'Arterial','LineWidth',2); hold on;
plot(t_4, A_4(:,2), 'DisplayName', 'Venous','LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (mg/L)');
title('Plasma Concentration');
set(gca,'FontSize',25);
legend show;
grid on;

% arterial and venous concentrations
subplot(2,2,2);
plot(t_4, A_4(:,3), 'DisplayName', 'Lung Tissue','LineWidth',2); hold on;
xlabel('Time (h)'); ylabel('Concentration (mg/L)');
title('Drug Concentration in Lung Tissue');
set(gca,'FontSize',25);
legend show;
grid on;

subplot(2,2,3);
% oral and gut compartments
plot(t_4, A_4(:,7), 'DisplayName', 'Absorbed Drug (AD)','LineWidth',2);hold on;
plot(t_4, A_4(:,6), 'DisplayName', 'Gut','LineWidth',2);
xlabel('Time (h)'); ylabel('Amount of Drug (mg)');
title('Drug Amount in Dose Compartments');
set(gca,'FontSize',25);
legend show;
grid on;

sgtitle('Rifampicin - Mouse Day 1 - Model 4 Ramachandran','FontSize',30);



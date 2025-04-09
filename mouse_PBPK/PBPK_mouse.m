
 
BW = 0.025; % bodyweight of mouse in kg
phys = loadPhysiologyMouse(BW);

%Rifampicin
fup = 0.029;    % Fraction unbound in plasma (Supp Table 1)
CLr = 19;       % Fraction of renal clearance fraction (Supp Table 1)
F = 1;          % Drug bioavailability (Supp Table 1)
ka = 0.3713;    % Rate of absorption (1/h) (Supp Table 3)It says the 'best fitted parameter'
CL = 0.037;     % Total body clearance (L/h) (Supp Table 3) It says the 'best fitted parameter'

%Bedaquiline
% fup = 0.029;    % Fraction unbound in plasma (Supp Table 1)
% CLr = 19;       % Renal clearance fraction (Supp Table 1)
% F = 1;          % Drug bioavailability (Supp Table 1)
% ka = 0.3713;    % Rate of absorption (1/h) (Supp Table 3)It says the 'best fitted parameter'
% CL = 0.037;     % Total body clearance (L/h) (Supp Table 3) It says the 'best fitted parameter'


A0 = zeros(1,9); A0(7) = 10*BW; % drug is 10 mg/kg

tspan = [0 24]; % 24h simulation
[t, y] = ode15s(@(t, y) odePBPK_Mouse(t, y, fup, CLr,F, ka, CL, phys), tspan, A0);

A=y;

% Plot various drug concentrations
figure;

subplot(2,2,1);

% oral and gut compartments
plot(t, A(:,7), 'DisplayName', 'Absorbed Drug (AD)','LineWidth',2);hold on;
plot(t, A(:,6), 'DisplayName', 'Gut','LineWidth',2);
xlabel('Time (h)'); ylabel('Amount of Drug (mg)');
title('Drug Amount in Dose Compartments');
set(gca,'FontSize',25);
legend show;
grid on;

% lung compartments
subplot(2,2,2);
plot(t, A(:,3)/phys.V.Lu, 'DisplayName', 'Lung Tissue','LineWidth',2); hold on;
xlabel('Time (h)'); ylabel('Concentration (mg/L)');
title('Drug Concentration in Lung Tissue');
set(gca,'FontSize',25);
legend show;
grid on;

% arterial and venous concentrations
subplot(2,2,3);
plot(t, A(:,1) / phys.V.A, 'DisplayName', 'Arterial','LineWidth',2); hold on;
plot(t, A(:,2) / phys.V.V, 'DisplayName', 'Venous','LineWidth',2);
xlabel('Time (h)'); ylabel('Concentration (mg/L)');
title('Plasma Concentration');
set(gca,'FontSize',25);
legend show;
grid on;

sgtitle('Rifampicin','FontSize',30);


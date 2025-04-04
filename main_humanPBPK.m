
% Load parameters
params = parameters_humanPBPK();

% Initial conditions
y0 = zeros(1, 18); % Adjust the size based on the number of compartments
y0(18) = (600); % Initial dose in the gut lumen of RIF (mg) **this should probably be changed to ug

% Time span for simulation
tStart = 0;
tEnd = 24;
% tspan = [tStart tEnd];
% n_timepoints = 48;
tspan = tStart:0.5:tEnd;
opts2 = odeset('NonNegative',1:18);

% Solve ODEs using ode15s
[t, y] = ode45(@(t, y) ODEs_human(t, y, params), tspan, y0,opts2);

% Plot results
figure;
plot(t, y(:, 18), '-o'); % Plot lung concentration as an example
xlabel('Time (hours)');
ylabel('Drug Concentration in Pleura (ug/mL)');
title('Drug Concentration in Pleura Over Time');

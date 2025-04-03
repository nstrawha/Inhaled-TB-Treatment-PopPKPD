
% Load parameters
params = parameters_humanPBPK();

% Initial conditions
y0 = zeros(1, 18); % Adjust the size based on the number of compartments
y0(18) = (600*1000); % Initial dose in the gut lumen of RIg (mg) conver to ug

% Time span for simulation
tStart = 0;
tEnd = 30;
tspan = [tStart tEnd];

% Solve ODEs using ode15s
[t, y] = ode15s(@(t, y) ODEs(t, y, params), tspan, y0);

% Plot results
figure;
plot(t, y(:, 14), '-o'); % Plot lung concentration as an example
xlabel('Time (hours)');
ylabel('Drug Concentration in Pleura (ug/mL)');
title('Drug Concentration in Pleura Over Time');

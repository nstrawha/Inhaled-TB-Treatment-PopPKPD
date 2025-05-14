%% This is the main script to execute the Lung Sim ODEs we were given by Denise Kirscher

%% Section 1: - read in parameter file and run LHS 
% Read in parameters from CSV file in this folder
ODEParamsRead = readtable('ODEParams.csv');

% Convert Min and Max to numeric array for easy indexing
ODEParamsReadArray = table2array(ODEParamsRead(:, 2:3));  % Columns 2 = Min, 3 = Max

%% LHS MATRIX  %%
runs=100;
numParams = size(ODEParamsReadArray, 1); % Number of parameters
LHSMatrix = zeros(runs, numParams);      % Preallocate

for i = 1:numParams
    xmin = ODEParamsReadArray(i, 1);
    xmax = ODEParamsReadArray(i, 2);
    xmed = median([xmin, xmax]);   % Really only necessary for normal distribution so just choosing the median for now
    LHSMatrix(:, i) = LHS_Call(xmin, xmed, xmax, 0, runs, 'unif');
end

%% Section 2: - execute ODEs for parameter sets in LHS matrix
t=0;
tEnd = 200; % days
tStart=t;
opts = odeset('NonNegative',1);
% opts is likely for non-negative (said by Christian Michael a post doc in the lab Feb 2025)

% Preallocate results
T_all = cell(runs,1); % to store time vectors
Y_all = cell(runs,1); % to store solutions

% start with 1 infected mac, 1 intracellular mtb and 1 activated mac
% can start with 0 everything else 
y0 = zeros(1,16); y0(2) = 1; y0(3) =1; y0(11)=1;

paramNames = ODEParamsRead.ParameterName; % Cell array of strings
for i = 1:runs
    % Assemble a struct with the parameter values for this run
    ODEParams = struct();
    for j = 1:length(paramNames)
        param = matlab.lang.makeValidName(paramNames{j});
        ODEParams.(param) = LHSMatrix(i, j);
    end

    % Solve ODE
    [tTemp, yTemp] = ode15s(@(t,y) ODEquations(t, y, ODEParams), [tStart tEnd], y0, opts);

    % Store result
    T_all{i} = tTemp;
    Y_all{i} = real(yTemp);
end

% Order of ODEs 
names= {'MR';... % resting macrophages
    'MI';...    % infected macrophages
    'MA';...    % activated macrophages
    'T0';...    % T cells
    'T1';...    % T cells
    'T2';...    % T cells
    'T80';...   % T cells
    'TC';...    % T cells
    'T8';...    % T cells
    'BE';...    % extracellular bacteria
    'BI';...    % intracellular bacteria
    'IG';...    % IG 
    'I12';...   % IL-12
    'I10';...   % IL-10
    'I4';...    % IL-4
    'TNF'};     % TNF


%% Plots

fig = figure();
fig.Position = [00 00 1920 1080];
set(0,'DefaultFigureWindowStyle','docked');
tiledlayout(4,4);
% Plot each variable
for j=1:16
    nexttile
    for i=1:runs
        plot(T_all{i}, Y_all{i}(:, j),'LineWidth',2);
        title(names{j});
        hold on;
    end
end

% Recreate Figure 4 from Wessler et al 2020. Assuming these are the
% variables that add up to Total CFU, Total Macrophage and Total CD3s
for i=1:runs
    CFU{i} = Y_all{i}(:,10) + Y_all{i}(:,11);
    Macs{i} = Y_all{i}(:,1) + Y_all{i}(:,2)+ Y_all{i}(:,3);
    CD3Tcells{i}= Y_all{i}(:,4) + Y_all{i}(:,5)+ Y_all{i}(:,6) + Y_all{i}(:,7) + Y_all{i}(:,8)+ Y_all{i}(:,9);
end

fig = figure();
fig.Position = [00 00 1920 1080];
set(0,'DefaultFigureWindowStyle','docked');
tiledlayout(1,3)

% Fig 4A - Total CFU
nexttile
for i=1:runs
    plot(T_all{i}, CFU{i},'LineWidth',2);
    set(gca, 'YScale', 'log');
    title('Total CFU');
    hold on;
end

% Fig 4B -  Total Macrophages
nexttile
for i=1:runs
    plot(T_all{i}, Macs{i},'LineWidth',2);
    set(gca, 'YScale', 'log');
    title('Total Macrophages');
    hold on;
end

% Fig 4C -  Total CD3 T cells
nexttile
for i=1:runs
    plot(T_all{i}, CD3Tcells{i},'LineWidth',2);
    set(gca, 'YScale', 'log');
    title('Total CD3 T cells');
    hold on;
end

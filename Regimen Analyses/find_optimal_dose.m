%% Find Optimal Dose
% Script adapted from run_PBPK to find optimal dose for RIF
% June 2025

clc;
clearvars;

addpath("Oral_Dose_ODEs/");
addpath("Lung_Dose_ODEs/");
addpath("Methods/");
addpath("Methods/PK_Metrics/");
addpath("Methods/Parameters/");

%% Set parameters

% Patient parameters
body_weight = 70; % kg
phys = loadPhysiology(body_weight);
kF_RIF = 0.252; % gut transit rate

% RIF parameters
tstep_RIF = 0.01;

ka_oral_RIF = 1.08;     % absorption rate [1/h]
kdiss_lung_RIF = 50;    % dissolution rate [1/h] from Himstedt et al.
CL_RIF = 7.86;          % systemic clearance [L/h]
fR_RIF = 0.1830;        % fractional renal clearance
kr_RIF = 0.17;          % gut reabsorption rate [1/h]
kmuc_RIF = 0;           % mucociliary clearance [1/h] from Himstedt et al. (not implemented)
br_frac_RIF = 9/49;     % proportion of RIF absorbed by bronchi from Himstedt et al.
effRB_RIF = 2.56;       % bronchi efflux ratio from Himstedt et al.
effRA_RIF = 3.67;       % alveolar efflux ratio from Himstedt et al.
pt_RIF = loadPartitionCoefficients("RIF");

% model compts      name                index in concentration array
compt_list_RIF =   ["Plasma";           % 1
                    "Arterial Blood";   % 2
                    "Lung";             % 3
                    "Pleura";           % 4
                    "Brain";            % 5
                    "Adipose Tissue";   % 6
                    "Heart";            % 7
                    "Muscle";           % 8
                    "Skin";             % 9
                    "Other Tissue";     % 10
                    "Bone";             % 11
                    "Spleen";           % 12
                    "Kidney";           % 13
                    "Gut";              % 14
                    "Liver";            % 15
                    "Lymph Node";       % 16
                    "Gut Lumen";        % 17
                    "ELFb"              % 18
                    "ELFa"              % 19
                    "Absorption"];      % 20;

ncompts_total_RIF = length(compt_list_RIF);


%% Set up plotting and calculations

ndays_RIF = 8; % day to calculate ODEs metrics for

% Decide which compartments to calculate
compts_to_plot_RIF = {"Lung", "Plasma", "Pleura", "Lymph Node", "Liver", "Kidney"};

%% Solve ODEs and record data
% set dosing info
min_lungdose_RIF = 100; % mg
max_lungdose_RIF = 1000; % mg
lungdose_incr_RIF = 50; % mg/increment
min_lungdose_freq_RIF = 1; % /day
max_lungdose_freq_RIF = 4; % /day

oral_dose_RIF = 600;    % mg
oral_dose_freq_RIF = 1; % 1x daily

% initialize cell arrays to record data
data_format = cell((max_lungdose_RIF - min_lungdose_RIF) / lungdose_incr_RIF + 3, (max_lungdose_freq_RIF - min_lungdose_freq_RIF) + 3);
doses_RIF = min_lungdose_RIF:lungdose_incr_RIF:max_lungdose_RIF;
freqs_RIF = min_lungdose_freq_RIF:max_lungdose_freq_RIF;

for row_idx = 1:length(doses_RIF)
    data_format{2 + row_idx, 2} = doses_RIF(row_idx);
end

for col_idx = 1:length(freqs_RIF)
    data_format{2, 2 + col_idx} = freqs_RIF(col_idx);
end

% labels
data_format{1, 1} = "C_max; % different from oral dose";
data_format{1, 3} = "Dose (mg)";
data_format{3, 1} = "Frequency (doses/day)";

% set up cells for all compartments
table_cells = cell(length(compts_to_plot_RIF), 1);
for compt_idx = 1:length(compts_to_plot_RIF)
    table_cells{compt_idx} = data_format;
end

%% iterate through dosing regimens
for freq_idx = min_lungdose_freq_RIF:max_lungdose_freq_RIF
    lung_dose_freq_RIF = freq_idx;

    for dose_idx = min_lungdose_RIF:lungdose_incr_RIF:max_lungdose_RIF
        lung_dose_RIF = dose_idx;

        % package params                % index
        params_RIF =   {oral_dose_RIF;  % 1
                        lung_dose_RIF;  % 2
                        ka_oral_RIF;    % 3
                        CL_RIF;         % 4
                        fR_RIF;         % 5
                        kr_RIF;         % 6
                        pt_RIF;         % 7
                        phys;           % 8
                        kF_RIF;         % 9
                        kdiss_lung_RIF; % 10
                        effRB_RIF;      % 11
                        effRA_RIF;      % 12
                        kmuc_RIF;       % 13
                        br_frac_RIF;    % 14
                        tstep_RIF};     % 15

        [t_oraldose_RIF, C_oraldose_RIF, t_lungdose_RIF, C_lungdose_RIF] = solveODEs("RIF", params_RIF, ncompts_total_RIF, ndays_RIF, ...
                                                                                        oral_dose_freq_RIF, lung_dose_freq_RIF);
        % calculate AUC for requested compartments
        for compt_idx = 1:length(compts_to_plot_RIF)
            current_compt = compts_to_plot_RIF{compt_idx};
            [idx_to_calc, ~] = find(string(compt_list_RIF) == current_compt);

            starting_idx = (ndays_RIF * 24 / tstep_RIF) - 24 / tstep_RIF + 1;
            Cmax_lung = max(C_lungdose_RIF(starting_idx:end, idx_to_calc));
            Cmax_oral = max(C_oraldose_RIF(starting_idx:end, idx_to_calc));
            perc_diff = (Cmax_lung - Cmax_oral) / Cmax_oral * 100;
            perc_diff = round(perc_diff, 2);

            cell_for_compt = table_cells{compt_idx};
            cell_for_compt{dose_idx / lungdose_incr_RIF + 1, freq_idx + 2} = perc_diff;
            table_cells{compt_idx} = cell_for_compt;

        end
    end
end

%% write output sheets
for compt_idx = 1:length(table_cells)
    current_compt = table_cells{compt_idx};
    compt_name = compts_to_plot_RIF{compt_idx};

    writecell(current_compt, "Cmax_comparisons.xlsx", "Sheet", compt_name);
end
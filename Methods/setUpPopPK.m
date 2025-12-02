function setUpPopPK(drug, relevant_compts, toxic_compts, MICs_TB, AUC_target, Cmax_target, n_pts, n_days, oral_dose, lung_dose, oral_dose_freq, lung_dose_freq)
% SETUPPOPPK
% This function uses the equations and parameters found in Ramachandran &
% Gadgil, 2023, in order to compare an oral dose of TB medication to an
% inhaled/lung dose of TB medication. 
%
% PLOTS GENERATED:
% - Concentration-time courses of all patients in relevant compts.
% - 10th, 50th, and 90th percentiles of patient timecourses in relevant
%   compts.
% - MIC distribution and AUC/MIC probability of target attainment (PTA) 
%   comparison
% - MIC distribution and C_max/MIC PTA comparison
%
% FILES GENERATED:
% - popPK_analysis.xlsx: compares PK metrics across patient populations for
%   oral vs. lung dosing
% - popPK_CFRs.xlsx: contains estimated cumulative fractions of response 
%   (CFRs) for AUC/MIC and C_max/MIC PTA target values
%

% set wd to parent dir 
script_path = mfilename('fullpath'); 
script_dir = fileparts(script_path); 
parent_dir = fileparts(script_dir); 
cd(parent_dir);


%% Set parameters

% fixed model parameters
days_to_plot = 1;
tstep = 0.01;

br_frac = 9/49;  % proportion of drug absorbed by bronchi from Himstedt et al.

% model compts      name                index in concentration array
compt_list =       ["Plasma";           % 1
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
                    "ELFb";             % 18
                    "ELFa";             % 19
                    "Absorption"];      % 20;

ncompts_total = length(compt_list);

% initialize cell arrays to store params and output
Cs_oral_store = cell(1, n_pts);
Cs_lung_store = cell(1, n_pts);


%% Iterate through patient equations

% sample physiological parameters
[bw_PD, vol_PDs, vol_frac_PDs, ...
    qc_PD, flow_PDs, flow_frac_PDs, ...
    params_store] = getParamPDs(drug, 0.2, 0.3); % CVs from Lyons et al.

% initialize parameter storage
vol_params_store = params_store{1};
flow_params_store = params_store{2};

ts = 0:tstep:(24 * n_days - tstep);

parfor pt_idx = 1:n_pts
    
    all_params = loadPhysParams(drug, bw_PD, vol_PDs, vol_frac_PDs, ...
                                         qc_PD, flow_PDs, flow_frac_PDs);

    % solve ODEs
    [C_oraldose_pt, ...
        C_lungdose_pt] = solveODEs(all_params, tstep, br_frac, ...
                                        ncompts_total, n_days, ...
                                        oral_dose, lung_dose, ...
                                        oral_dose_freq, lung_dose_freq);

    % store results
    Cs_oral_store{pt_idx} = C_oraldose_pt;
    Cs_lung_store{pt_idx} = C_lungdose_pt;

    % store parameters
    vol_params = all_params{1};
    flow_params = all_params{2};
    vol_params_store = [vol_params_store; vol_params];
    flow_params_store = [flow_params_store; flow_params];

    % track progress
    disp(append("pt. #", num2str(pt_idx), " calculated..."))

end

% repackage parameter storage
params_store{1} = vol_params_store;
params_store{2} = flow_params_store;


%% Plot resulting timecourses

% iterate through compartments
for compt_idx = 1:length(relevant_compts)
    compt = relevant_compts{compt_idx};
    [idx_to_plot, ~] = find(string(compt_list) == compt);

    % pull all concentration timecourses for the current compartment
    current_Cs_orals = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_oral_store, "UniformOutput", false));
    current_Cs_lungs = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_lung_store, "UniformOutput", false));
    size(current_Cs_orals)
    plotTimecourses(drug, oral_dose, lung_dose, ...
                    compt, n_days, days_to_plot, ...
                    oral_dose_freq, lung_dose_freq, ...
                    ts, current_Cs_orals, current_Cs_lungs);
end


%% Analysis

% compare PK metrics
[AUCs_oral, ...
    AUCs_lung, ...
    Cmaxs_oral, ...
    Cmaxs_lung] = trackPKMetrics(drug, compt_list, relevant_compts, toxic_compts, ...
                                        oral_dose, oral_dose_freq, ...
                                        lung_dose, lung_dose_freq, ...
                                        n_days, tstep, ts, n_pts, ...
                                        Cs_oral_store, Cs_lung_store);


%% Calculate PTAs for oral and lung dosing

% record nontoxic compartments to calculate PTA for
compts_to_calc = "Lung";

if AUC_target ~= 0
    % calculate and compare PTAs
    plotPTAs(drug, "AUC/MIC", compts_to_calc, ...
                AUCs_oral, AUCs_lung, ...
                Cmaxs_oral, Cmaxs_lung, ...
                MICs_TB, AUC_target, ...
                oral_dose, oral_dose_freq, ...
                lung_dose, lung_dose_freq, ...
                n_pts, n_days);
end

if Cmax_target ~= 0
    plotPTAs(drug, "Cmax/MIC", compts_to_calc, ...
                AUCs_oral, AUCs_lung, ...
                Cmaxs_oral, Cmaxs_lung, ...
                MICs_TB, Cmax_target, ...
                oral_dose, oral_dose_freq, ...
                lung_dose, lung_dose_freq, ...
                n_pts, n_days);
end
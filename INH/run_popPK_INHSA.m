%% Run PopPK (INHSA)
% 
% DESCRIPTION:
% This script uses the equations and parameters found in Ramachandran &
% Gadgil, 2023, in order to compare an oral dose of isoniazid (slow 
% acetylator) to an inhaled/lung dose of isoniazid. 
% 
% The dosing regimen of each may be altered with oral/lung_dose_INHSA and
% oral/lung_dose_freq_INHSA. The compartments which will be plotted and for
% which PK metrics (AUC_24, C_avg, and C_max) will be calculated may be
% altered with relevant_compts_INHSA.
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
% AUTHORS:
% Noah Strawhacker (nstrawha@purdue.edu)
% Alexis Hoerter

clc;
clearvars;

% set wd to parent dir 
script_path = mfilename('fullpath'); 
script_dir = fileparts(script_path); 
parent_dir = fileparts(script_dir); 
cd(parent_dir);

addpath("INHSA/");
addpath("Methods/");


%% Set parameters

% fixed model parameters
n_pts_INHSA = 1;
n_days_INHSA = 8; % > TODO is steady state
days_to_plot_INHSA = 1;
relevant_compts_INHSA = {"Plasma", "Pleura", "Lung"};

oral_dose_INHSA = 300;    % mg
lung_dose_INHSA = 300;    % mg
oral_dose_freq_INHSA = 1; % doses/day
lung_dose_freq_INHSA = 1; % doses/day

tstep_INHSA = 0.01;

effRB_INHSA = 2.56;    % bronchi efflux ratio from Himstedt et al.
effRA_INHSA = 3.67;    % alveolar efflux ratio from Himstedt et al.
br_frac_INHSA = 9/49;  % proportion of INHSA absorbed by bronchi from Himstedt et al.

% model compts      name                index in concentration array
compt_list_INHSA =   ["Plasma";           % 1
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

ncompts_total_INHSA = length(compt_list_INHSA);
toxic_compts_INHSA = ["Liver", "Kidney"];

% initialize cell arrays to store params and output
Cs_oral_store_INHSA = cell(1, n_pts_INHSA);
Cs_lung_store_INHSA = cell(1, n_pts_INHSA);


%% Iterate through patient equations

% sample physiological parameters
[bw_PD, vol_PDs_INHSA, vol_frac_PDs_INHSA, ...
    qc_PD, flow_PDs_INHSA, flow_frac_PDs_INHSA, ...
    params_store_INHSA] = getParamPDs("INHSA", 0.00000000000001, 0.0000000000001); % CVs from Lyons et al.

% initialize parameter storage
vol_params_store_INHSA = params_store_INHSA{1};
flow_params_store_INHSA = params_store_INHSA{2};

ts_INHSA = 0:tstep_INHSA:(24 * n_days_INHSA - tstep_INHSA);

parfor pt_idx = 1:n_pts_INHSA
    
    all_params_INHSA = loadPhysParams("INHSA", bw_PD, vol_PDs_INHSA, vol_frac_PDs_INHSA, ...
                                         qc_PD, flow_PDs_INHSA, flow_frac_PDs_INHSA);

    % solve ODEs
    [C_oraldose_INHSA_pt, ...
        C_lungdose_INHSA_pt] = solveODEs(all_params_INHSA, tstep_INHSA, br_frac_INHSA, ...
                                        effRB_INHSA, effRA_INHSA, ...
                                        ncompts_total_INHSA, n_days_INHSA, ...
                                        oral_dose_INHSA, lung_dose_INHSA, ...
                                        oral_dose_freq_INHSA, lung_dose_freq_INHSA);

    % store results
    Cs_oral_store_INHSA{pt_idx} = C_oraldose_INHSA_pt;
    Cs_lung_store_INHSA{pt_idx} = C_lungdose_INHSA_pt;

    % store parameters
    vol_params_INHSA = all_params_INHSA{1};
    flow_params_INHSA = all_params_INHSA{2};
    vol_params_store_INHSA = [vol_params_store_INHSA; vol_params_INHSA];
    flow_params_store_INHSA = [flow_params_store_INHSA; flow_params_INHSA];

    % track progress
    disp(append("pt. #", num2str(pt_idx), " calculated..."))

end

% repackage parameter storage
params_store_INHSA{1} = vol_params_store_INHSA;
params_store_INHSA{2} = flow_params_store_INHSA;


%% Plot resulting timecourses

% iterate through compartments
for compt_idx = 1:length(relevant_compts_INHSA)
    compt = relevant_compts_INHSA{compt_idx};
    [idx_to_plot, ~] = find(string(compt_list_INHSA) == compt);

    % pull all concentration timecourses for the current compartment
    current_Cs_orals = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_oral_store_INHSA, "UniformOutput", false));
    current_Cs_lungs = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_lung_store_INHSA, "UniformOutput", false));

    plotTimecourses("INHSA", oral_dose_INHSA, lung_dose_INHSA, ...
                    compt, n_days_INHSA, days_to_plot_INHSA, ...
                    oral_dose_freq_INHSA, lung_dose_freq_INHSA, ...
                    ts_INHSA, current_Cs_orals, current_Cs_lungs);
end


%% Analysis

% compare PK metrics
[AUCs_oral_INHSA, ...
    AUCs_lung_INHSA, ...
    Cmaxs_oral_INHSA, ...
    Cmaxs_lung_INHSA] = trackPKMetrics("INHSA", compt_list_INHSA, relevant_compts_INHSA, toxic_compts_INHSA, ...
                                        oral_dose_INHSA, oral_dose_freq_INHSA, ...
                                        lung_dose_INHSA, lung_dose_freq_INHSA, ...
                                        n_days_INHSA, tstep_INHSA, ts_INHSA, n_pts_INHSA, ...
                                        Cs_oral_store_INHSA, Cs_lung_store_INHSA);


%% Calculate PTAs for oral and lung dosing

% MIC dist formatted as [conc., # isolates]
MICs_TB_INHSA =  {[0.031, 28];
                [0.062, 47];
                [0.120, 51];
                [0.250, 76];
                [0.500, 81];
                [1.000, 68];
                [2.000, 0];
                [4.000, 5];
                [8.000, 0];
                [16.00, 0]};

INHFA_AUC_target = 567; % AUC/MIC, Alffenaar et al.

% record nontoxic compartments to calculate PTA for
nontoxic_compts_INHSA = [];
for compt_idx = 1:length(relevant_compts_INHSA)
    current_compt = relevant_compts_INHSA{compt_idx};

    if ~ismember(current_compt, toxic_compts_INHSA)
        nontoxic_compts_INHSA = [nontoxic_compts_INHSA, current_compt];
    end

end

% calculate and compare PTAs
plotPTAs("INHSA", "AUC/MIC", nontoxic_compts_INHSA, ...
            AUCs_oral_INHSA, AUCs_lung_INHSA, ...
            Cmaxs_oral_INHSA, Cmaxs_lung_INHSA, ...
            MICs_TB_INHSA, INHSA_AUC_target, ...
            oral_dose_INHSA, oral_dose_freq_INHSA, ...
            lung_dose_INHSA, lung_dose_freq_INHSA, ...
            n_pts_INHSA, n_days_INHSA);
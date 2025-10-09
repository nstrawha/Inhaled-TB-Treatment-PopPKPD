%% Run PopPK (INHFA)
% 
% DESCRIPTION:
% This script uses the equations and parameters found in Ramachandran &
% Gadgil, 2023, in order to compare an oral dose of isoniazid (fast 
% acetylator) to an inhaled/lung dose of isoniazid. 
% 
% The dosing regimen of each may be altered with oral/lung_dose_INHFA and
% oral/lung_dose_freq_INHFA. The compartments which will be plotted and for
% which PK metrics (AUC_24, C_avg, and C_max) will be calculated may be
% altered with relevant_compts_INHFA.
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

addpath("INHFA/");
addpath("Methods/");


%% Set parameters

% fixed model parameters
n_pts_INHFA = 1;
n_days_INHFA = 8; % > TODO is steady state
days_to_plot_INHFA = 1;
relevant_compts_INHFA = {"Plasma", "Pleura", "Lung"};

oral_dose_INHFA = 300;    % mg
lung_dose_INHFA = 300;    % mg
oral_dose_freq_INHFA = 1; % doses/day
lung_dose_freq_INHFA = 1; % doses/day

tstep_INHFA = 0.01;

effRB_INHFA = 2.56;    % bronchi efflux ratio from Himstedt et al.
effRA_INHFA = 3.67;    % alveolar efflux ratio from Himstedt et al.
br_frac_INHFA = 9/49;  % proportion of INHFA absorbed by bronchi from Himstedt et al.

% model compts      name                index in concentration array
compt_list_INHFA =   ["Plasma";           % 1
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

ncompts_total_INHFA = length(compt_list_INHFA);
toxic_compts_INHFA = ["Liver", "Kidney"];

% initialize cell arrays to store params and output
Cs_oral_store_INHFA = cell(1, n_pts_INHFA);
Cs_lung_store_INHFA = cell(1, n_pts_INHFA);


%% Iterate through patient equations

% sample physiological parameters
[bw_PD, vol_PDs_INHFA, vol_frac_PDs_INHFA, ...
    qc_PD, flow_PDs_INHFA, flow_frac_PDs_INHFA, ...
    params_store_INHFA] = getParamPDs("INHFA", 0.00000000000001, 0.0000000000001); % CVs from Lyons et al.

% initialize parameter storage
vol_params_store_INHFA = params_store_INHFA{1};
flow_params_store_INHFA = params_store_INHFA{2};

ts_INHFA = 0:tstep_INHFA:(24 * n_days_INHFA - tstep_INHFA);

parfor pt_idx = 1:n_pts_INHFA
    
    all_params_INHFA = loadPhysParams("INHFA", bw_PD, vol_PDs_INHFA, vol_frac_PDs_INHFA, ...
                                         qc_PD, flow_PDs_INHFA, flow_frac_PDs_INHFA);

    % solve ODEs
    [C_oraldose_INHFA_pt, ...
        C_lungdose_INHFA_pt] = solveODEs(all_params_INHFA, tstep_INHFA, br_frac_INHFA, ...
                                        effRB_INHFA, effRA_INHFA, ...
                                        ncompts_total_INHFA, n_days_INHFA, ...
                                        oral_dose_INHFA, lung_dose_INHFA, ...
                                        oral_dose_freq_INHFA, lung_dose_freq_INHFA);

    % store results
    Cs_oral_store_INHFA{pt_idx} = C_oraldose_INHFA_pt;
    Cs_lung_store_INHFA{pt_idx} = C_lungdose_INHFA_pt;

    % store parameters
    vol_params_INHFA = all_params_INHFA{1};
    flow_params_INHFA = all_params_INHFA{2};
    vol_params_store_INHFA = [vol_params_store_INHFA; vol_params_INHFA];
    flow_params_store_INHFA = [flow_params_store_INHFA; flow_params_INHFA];

    % track progress
    disp(append("pt. #", num2str(pt_idx), " calculated..."))

end

% repackage parameter storage
params_store_INHFA{1} = vol_params_store_INHFA;
params_store_INHFA{2} = flow_params_store_INHFA;


%% Plot resulting timecourses

% iterate through compartments
for compt_idx = 1:length(relevant_compts_INHFA)
    compt = relevant_compts_INHFA{compt_idx};
    [idx_to_plot, ~] = find(string(compt_list_INHFA) == compt);

    % pull all concentration timecourses for the current compartment
    current_Cs_orals = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_oral_store_INHFA, "UniformOutput", false));
    current_Cs_lungs = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_lung_store_INHFA, "UniformOutput", false));

    plotTimecourses("INHFA", oral_dose_INHFA, lung_dose_INHFA, ...
                    compt, n_days_INHFA, days_to_plot_INHFA, ...
                    oral_dose_freq_INHFA, lung_dose_freq_INHFA, ...
                    ts_INHFA, current_Cs_orals, current_Cs_lungs);
end


%% Analysis

% compare PK metrics
[AUCs_oral_INHFA, ...
    AUCs_lung_INHFA, ...
    Cmaxs_oral_INHFA, ...
    Cmaxs_lung_INHFA] = trackPKMetrics("INHFA", compt_list_INHFA, relevant_compts_INHFA, toxic_compts_INHFA, ...
                                        oral_dose_INHFA, oral_dose_freq_INHFA, ...
                                        lung_dose_INHFA, lung_dose_freq_INHFA, ...
                                        n_days_INHFA, tstep_INHFA, ts_INHFA, n_pts_INHFA, ...
                                        Cs_oral_store_INHFA, Cs_lung_store_INHFA);


%% Calculate PTAs for oral and lung dosing

% MIC dist formatted as [conc., # isolates]
MICs_TB_INHFA =  {[0.031, 28];
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
nontoxic_compts_INHFA = [];
for compt_idx = 1:length(relevant_compts_INHFA)
    current_compt = relevant_compts_INHFA{compt_idx};

    if ~ismember(current_compt, toxic_compts_INHFA)
        nontoxic_compts_INHFA = [nontoxic_compts_INHFA, current_compt];
    end

end

% calculate and compare PTAs
plotPTAs("INHFA", "AUC/MIC", nontoxic_compts_INHFA, ...
            AUCs_oral_INHFA, AUCs_lung_INHFA, ...
            Cmaxs_oral_INHFA, Cmaxs_lung_INHFA, ...
            MICs_TB_INHFA, INHFA_AUC_target, ...
            oral_dose_INHFA, oral_dose_freq_INHFA, ...
            lung_dose_INHFA, lung_dose_freq_INHFA, ...
            n_pts_INHFA, n_days_INHFA);
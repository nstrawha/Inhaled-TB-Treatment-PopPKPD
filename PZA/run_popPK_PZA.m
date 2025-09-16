%% Run PopPK (PZA)
% 
% DESCRIPTION:
% This script uses the equations and parameters found in Ramachandran &
% Gadgil, 2023, in order to compare an oral dose of pyrazinamide to an
% inhaled/lung dose of pyrazinamide. 
% 
% The dosing regimen of each may be altered with oral/lung_dose_PZA and
% oral/lung_dose_freq_PZA. The compartments which will be plotted and for
% which PK metrics (AUC_24, C_avg, and C_max) will be calculated may be
% altered with relevant_compts_PZA.
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

addpath("PZA/");
addpath("Methods/");


%% Set parameters

% fixed model parameters
n_pts_PZA = 1;
n_days_PZA = 1; % > TODO is steady state
days_to_plot_PZA = 1;
relevant_compts_PZA = {"Plasma"};

oral_dose_PZA = 2000;    % mg
lung_dose_PZA = 2000;    % mg
oral_dose_freq_PZA = 1;  % doses/day
lung_dose_freq_PZA = 1;  % doses/day

tstep_PZA = 0.01;

effRB_PZA = 2.56;    % bronchi efflux ratio from Himstedt et al.
effRA_PZA = 3.67;    % alveolar efflux ratio from Himstedt et al.
br_frac_PZA = 9/49;  % proportion of PZA absorbed by bronchi from Himstedt et al.

% model compts      name                index in concentration array
compt_list_PZA =   ["Plasma";           % 1
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

ncompts_total_PZA = length(compt_list_PZA);
toxic_compts_PZA = ["Liver", "Kidney"];

% initialize cell arrays to store params and output
Cs_oral_store_PZA = cell(1, n_pts_PZA);
Cs_lung_store_PZA = cell(1, n_pts_PZA);


%% Iterate through patient equations

% sample physiological parameters
[bw_PD, vol_PDs_PZA, vol_frac_PDs_PZA, ...
    qc_PD, flow_PDs_PZA, flow_frac_PDs_PZA, ...
    params_store_PZA] = getParamPDs("PZA", 0.00000000000001, 0.0000000000001); % CVs from Lyons et al.

% initialize parameter storage
vol_params_store_PZA = params_store_PZA{1};
flow_params_store_PZA = params_store_PZA{2};

ts_PZA = 0:tstep_PZA:(24 * n_days_PZA - tstep_PZA);

parfor pt_idx = 1:n_pts_PZA
    
    all_params_PZA = loadPhysParams("PZA", bw_PD, vol_PDs_PZA, vol_frac_PDs_PZA, ...
                                         qc_PD, flow_PDs_PZA, flow_frac_PDs_PZA);

    % solve ODEs
    [C_oraldose_PZA_pt, ...
        C_lungdose_PZA_pt] = solveODEs(all_params_PZA, tstep_PZA, br_frac_PZA, ...
                                        effRB_PZA, effRA_PZA, ...
                                        ncompts_total_PZA, n_days_PZA, ...
                                        oral_dose_PZA, lung_dose_PZA, ...
                                        oral_dose_freq_PZA, lung_dose_freq_PZA);

    % store results
    Cs_oral_store_PZA{pt_idx} = C_oraldose_PZA_pt;
    Cs_lung_store_PZA{pt_idx} = C_lungdose_PZA_pt;

    % store parameters
    vol_params_PZA = all_params_PZA{1};
    flow_params_PZA = all_params_PZA{2};
    vol_params_store_PZA = [vol_params_store_PZA; vol_params_PZA];
    flow_params_store_PZA = [flow_params_store_PZA; flow_params_PZA];

    % track progress
    disp(append("pt. #", num2str(pt_idx), " calculated..."))

end

% repackage parameter storage
params_store_PZA{1} = vol_params_store_PZA;
params_store_PZA{2} = flow_params_store_PZA;


%% Plot resulting timecourses

% iterate through compartments
for compt_idx = 1:length(relevant_compts_PZA)
    compt = relevant_compts_PZA{compt_idx};
    [idx_to_plot, ~] = find(string(compt_list_PZA) == compt);

    % pull all concentration timecourses for the current compartment
    current_Cs_orals = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_oral_store_PZA, "UniformOutput", false));
    current_Cs_lungs = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_lung_store_PZA, "UniformOutput", false));

    plotTimecourses("PZA", oral_dose_PZA, lung_dose_PZA, ...
                    compt, n_days_PZA, days_to_plot_PZA, ...
                    oral_dose_freq_PZA, lung_dose_freq_PZA, ...
                    ts_PZA, current_Cs_orals, current_Cs_lungs);
end


%% Analysis

% compare PK metrics
[AUCs_oral_PZA, ...
    AUCs_lung_PZA, ...
    Cmaxs_oral_PZA, ...
    Cmaxs_lung_PZA] = trackPKMetrics("PZA", compt_list_PZA, relevant_compts_PZA, toxic_compts_PZA, ...
                                        oral_dose_PZA, oral_dose_freq_PZA, ...
                                        lung_dose_PZA, lung_dose_freq_PZA, ...
                                        n_days_PZA, tstep_PZA, ts_PZA, n_pts_PZA, ...
                                        Cs_oral_store_PZA, Cs_lung_store_PZA);


%% Calculate PTAs for oral and lung dosing

% MIC dist formatted as [conc., # isolates]
MICs_TB_PZA =  {[0.031, 28];
                [0.062, 47];
                [0.120, 51];
                [0.250, 76];
                [0.500, 81];
                [1.000, 68];
                [2.000, 0];
                [4.000, 5];
                [8.000, 0];
                [16.00, 0]};

PZA_AUC_target = 120; % AUC/MIC, from Gumbo et al.

% record nontoxic compartments to calculate PTA for
nontoxic_compts_PZA = [];
for compt_idx = 1:length(relevant_compts_PZA)
    current_compt = relevant_compts_PZA{compt_idx};

    if ~ismember(current_compt, toxic_compts_PZA)
        nontoxic_compts_PZA = [nontoxic_compts_PZA, current_compt];
    end

end

% calculate and compare PTAs
plotPTAs("PZA", "AUC/MIC", nontoxic_compts_PZA, ...
            AUCs_oral_PZA, AUCs_lung_PZA, ...
            Cmaxs_oral_PZA, Cmaxs_lung_PZA, ...
            MICs_TB_PZA, PZA_AUC_target, ...
            oral_dose_PZA, oral_dose_freq_PZA, ...
            lung_dose_PZA, lung_dose_freq_PZA, ...
            n_pts_PZA, n_days_PZA);
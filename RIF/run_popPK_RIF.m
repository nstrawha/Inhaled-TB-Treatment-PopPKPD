%% Run PopPK (RIF)
% 
% DESCRIPTION:
% This script uses the equations and parameters found in Ramachandran &
% Gadgil, 2023, in order to compare an oral dose of rifampin to an
% inhaled/lung dose of rifampin. 
% 
% The dosing regimen of each may be altered with oral/lung_dose_RIF and
% oral/lung_dose_freq_RIF. The compartments which will be plotted and for
% which PK metrics (AUC_24, C_avg, and C_max) will be calculated may be
% altered with relevant_compts_RIF.
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

addpath("RIF/");
addpath("Methods/");


%% Set parameters

% fixed model parameters
n_pts_RIF = 1;
n_days_RIF = 4; % > 3 is steady state
days_to_plot_RIF = 1;
relevant_compts_RIF = {"Plasma"};

oral_dose_RIF = 600;    % mg
lung_dose_RIF = 600;    % mg
oral_dose_freq_RIF = 1; % doses/day
lung_dose_freq_RIF = 1; % doses/day

tstep_RIF = 0.01;

effRB_RIF = 2.56;    % bronchi efflux ratio from Himstedt et al.
effRA_RIF = 3.67;    % alveolar efflux ratio from Himstedt et al.
br_frac_RIF = 9/49;  % proportion of RIF absorbed by bronchi from Himstedt et al.

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
                    "ELFb";             % 18
                    "ELFa";             % 19
                    "Absorption"];      % 20;

ncompts_total_RIF = length(compt_list_RIF);
toxic_compts_RIF = ["Liver", "Kidney"];

% initialize cell arrays to store params and output
Cs_oral_store_RIF = cell(1, n_pts_RIF);
Cs_lung_store_RIF = cell(1, n_pts_RIF);


%% Iterate through patient equations

% sample physiological parameters
[bw_PD, vol_PDs_RIF, vol_frac_PDs_RIF, ...
    qc_PD, flow_PDs_RIF, flow_frac_PDs_RIF, ...
    params_store_RIF] = getParamPDs("RIF", 0.00000000000001, 0.0000000000001); % CVs from Lyons et al.

% initialize parameter storage
vol_params_store_RIF = params_store_RIF{1};
flow_params_store_RIF = params_store_RIF{2};

ts_RIF = 0:tstep_RIF:(24 * n_days_RIF - tstep_RIF);

parfor pt_idx = 1:n_pts_RIF
    
    all_params_RIF = loadPhysParams("RIF", bw_PD, vol_PDs_RIF, vol_frac_PDs_RIF, ...
                                         qc_PD, flow_PDs_RIF, flow_frac_PDs_RIF);

    % solve ODEs
    [C_oraldose_RIF_pt, ...
        C_lungdose_RIF_pt] = solveODEs(all_params_RIF, tstep_RIF, br_frac_RIF, ...
                                        effRB_RIF, effRA_RIF, ...
                                        ncompts_total_RIF, n_days_RIF, ...
                                        oral_dose_RIF, lung_dose_RIF, ...
                                        oral_dose_freq_RIF, lung_dose_freq_RIF);

    % store results
    Cs_oral_store_RIF{pt_idx} = C_oraldose_RIF_pt;
    Cs_lung_store_RIF{pt_idx} = C_lungdose_RIF_pt;

    % store parameters
    vol_params_RIF = all_params_RIF{1};
    flow_params_RIF = all_params_RIF{2};
    vol_params_store_RIF = [vol_params_store_RIF; vol_params_RIF];
    flow_params_store_RIF = [flow_params_store_RIF; flow_params_RIF];

    % track progress
    disp(append("pt. #", num2str(pt_idx), " calculated..."))

end

% repackage parameter storage
params_store_RIF{1} = vol_params_store_RIF;
params_store_RIF{2} = flow_params_store_RIF;


%% Plot resulting timecourses

% iterate through compartments
for compt_idx = 1:length(relevant_compts_RIF)
    compt = relevant_compts_RIF{compt_idx};
    [idx_to_plot, ~] = find(string(compt_list_RIF) == compt);

    % pull all concentration timecourses for the current compartment
    current_Cs_orals = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_oral_store_RIF, "UniformOutput", false));
    current_Cs_lungs = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_lung_store_RIF, "UniformOutput", false));

    plotTimecourses("RIF", oral_dose_RIF, lung_dose_RIF, ...
                    compt, n_days_RIF, days_to_plot_RIF, ...
                    oral_dose_freq_RIF, lung_dose_freq_RIF, ...
                    ts_RIF, current_Cs_orals, current_Cs_lungs);
end


%% Analysis

% compare PK metrics
[AUCs_oral_RIF, ...
    AUCs_lung_RIF, ...
    Cmaxs_oral_RIF, ...
    Cmaxs_lung_RIF] = trackPKMetrics("RIF", compt_list_RIF, relevant_compts_RIF, toxic_compts_RIF, ...
                                        oral_dose_RIF, oral_dose_freq_RIF, ...
                                        lung_dose_RIF, lung_dose_freq_RIF, ...
                                        n_days_RIF, tstep_RIF, ts_RIF, n_pts_RIF, ...
                                        Cs_oral_store_RIF, Cs_lung_store_RIF);


%% Calculate PTAs for oral and lung dosing

% MIC dist formatted as [conc., # isolates]
MICs_TB_RIF =  {[0.031, 28];
                [0.062, 47];
                [0.120, 51];
                [0.250, 76];
                [0.500, 81];
                [1.000, 68];
                [2.000, 0];
                [4.000, 5];
                [8.000, 0];
                [16.00, 0]};

RIF_AUC_target = 271; % AUC/MIC, from Jayaram et al.
RIF_Cmax_target = 175; % Cmax/MIC, from Gumbo et al.

% record nontoxic compartments to calculate PTA for
nontoxic_compts_RIF = [];
for compt_idx = 1:length(relevant_compts_RIF)
    current_compt = relevant_compts_RIF{compt_idx};

    if ~ismember(current_compt, toxic_compts_RIF)
        nontoxic_compts_RIF = [nontoxic_compts_RIF, current_compt];
    end

end

% calculate and compare PTAs
plotPTAs("RIF", nontoxic_compts_RIF, ...
            AUCs_oral_RIF, AUCs_lung_RIF, ...
            Cmaxs_oral_RIF, Cmaxs_lung_RIF, ...
            MICs_TB_RIF, ...
            RIF_AUC_target, RIF_Cmax_target, ...
            oral_dose_RIF, oral_dose_freq_RIF, ...
            lung_dose_RIF, lung_dose_freq_RIF, ...
            n_pts_RIF, n_days_RIF);
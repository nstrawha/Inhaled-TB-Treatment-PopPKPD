%% Run PopPK (EMB)
% 
% DESCRIPTION:
% This script uses the equations and parameters found in Ramachandran &
% Gadgil, 2023, in order to compare an oral dose of ethambutol to an
% inhaled/lung dose of ethambutol. 
% 
% The dosing regimen of each may be altered with oral/lung_dose_EMB and
% oral/lung_dose_freq_EMB. The compartments which will be plotted and for
% which PK metrics (AUC_24, C_avg, and C_max) will be calculated may be
% altered with relevant_compts_EMB.
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

addpath("EMB/");
addpath("Methods/");


%% Set parameters

% fixed model parameters
n_pts_EMB = 1;
n_days_EMB = 1; % > TODO is steady state
days_to_plot_EMB = 1;
relevant_compts_EMB = {"Lung"};

oral_dose_EMB = 1500;    % mg
lung_dose_EMB = 1500;    % mg
oral_dose_freq_EMB = 1;  % doses/day
lung_dose_freq_EMB = 1;  % doses/day

tstep_EMB = 0.01;

effRB_EMB = 2.56;    % bronchi efflux ratio from Himstedt et al.
effRA_EMB = 3.67;    % alveolar efflux ratio from Himstedt et al.
br_frac_EMB = 9/49;  % proportion of EMB absorbed by bronchi from Himstedt et al.

% model compts      name                index in concentration array
compt_list_EMB =   ["Plasma";           % 1
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

ncompts_total_EMB = length(compt_list_EMB);
toxic_compts_EMB = ["Liver", "Kidney"];

% initialize cell arrays to store params and output
Cs_oral_store_EMB = cell(1, n_pts_EMB);
Cs_lung_store_EMB = cell(1, n_pts_EMB);


%% Iterate through patient equations

% sample physiological parameters
[bw_PD, vol_PDs_EMB, vol_frac_PDs_EMB, ...
    qc_PD, flow_PDs_EMB, flow_frac_PDs_EMB, ...
    params_store_EMB] = getParamPDs("EMB", 0.00000000000001, 0.0000000000001); % CVs from Lyons et al.

% initialize parameter storage
vol_params_store_EMB = params_store_EMB{1};
flow_params_store_EMB = params_store_EMB{2};

ts_EMB = 0:tstep_EMB:(24 * n_days_EMB - tstep_EMB);

parfor pt_idx = 1:n_pts_EMB
    
    all_params_EMB = loadPhysParams("EMB", bw_PD, vol_PDs_EMB, vol_frac_PDs_EMB, ...
                                         qc_PD, flow_PDs_EMB, flow_frac_PDs_EMB);

    % solve ODEs
    [C_oraldose_EMB_pt, ...
        C_lungdose_EMB_pt] = solveODEs(all_params_EMB, tstep_EMB, br_frac_EMB, ...
                                        effRB_EMB, effRA_EMB, ...
                                        ncompts_total_EMB, n_days_EMB, ...
                                        oral_dose_EMB, lung_dose_EMB, ...
                                        oral_dose_freq_EMB, lung_dose_freq_EMB);

    % store results
    Cs_oral_store_EMB{pt_idx} = C_oraldose_EMB_pt;
    Cs_lung_store_EMB{pt_idx} = C_lungdose_EMB_pt;

    % store parameters
    vol_params_EMB = all_params_EMB{1};
    flow_params_EMB = all_params_EMB{2};
    vol_params_store_EMB = [vol_params_store_EMB; vol_params_EMB];
    flow_params_store_EMB = [flow_params_store_EMB; flow_params_EMB];

    % track progress
    disp(append("pt. #", num2str(pt_idx), " calculated..."))

end

% repackage parameter storage
params_store_EMB{1} = vol_params_store_EMB;
params_store_EMB{2} = flow_params_store_EMB;


%% Plot resulting timecourses

% iterate through compartments
for compt_idx = 1:length(relevant_compts_EMB)
    compt = relevant_compts_EMB{compt_idx};
    [idx_to_plot, ~] = find(string(compt_list_EMB) == compt);

    % pull all concentration timecourses for the current compartment
    current_Cs_orals = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_oral_store_EMB, "UniformOutput", false));
    current_Cs_lungs = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_lung_store_EMB, "UniformOutput", false));

    plotTimecourses("EMB", oral_dose_EMB, lung_dose_EMB, ...
                    compt, n_days_EMB, days_to_plot_EMB, ...
                    oral_dose_freq_EMB, lung_dose_freq_EMB, ...
                    ts_EMB, current_Cs_orals, current_Cs_lungs);
end


%% Analysis

% compare PK metrics
[AUCs_oral_EMB, ...
    AUCs_lung_EMB, ...
    Cmaxs_oral_EMB, ...
    Cmaxs_lung_EMB] = trackPKMetrics("EMB", compt_list_EMB, relevant_compts_EMB, toxic_compts_EMB, ...
                                        oral_dose_EMB, oral_dose_freq_EMB, ...
                                        lung_dose_EMB, lung_dose_freq_EMB, ...
                                        n_days_EMB, tstep_EMB, ts_EMB, n_pts_EMB, ...
                                        Cs_oral_store_EMB, Cs_lung_store_EMB);


%% Calculate PTAs for oral and lung dosing

% MIC dist formatted as [conc., # isolates]
MICs_TB_EMB =  {[0.031, 28];
                [0.062, 47];
                [0.120, 51];
                [0.250, 76];
                [0.500, 81];
                [1.000, 68];
                [2.000, 0];
                [4.000, 5];
                [8.000, 0];
                [16.00, 0]};

EMB_AUC_target = 271; % AUC/MIC, from Jayaram et al.
EMB_Cmax_target = 175; % Cmax/MIC, from Gumbo et al.

% record nontoxic compartments to calculate PTA for
nontoxic_compts_EMB = [];
for compt_idx = 1:length(relevant_compts_EMB)
    current_compt = relevant_compts_EMB{compt_idx};

    if ~ismember(current_compt, toxic_compts_EMB)
        nontoxic_compts_EMB = [nontoxic_compts_EMB, current_compt];
    end

end

% calculate and compare PTAs
plotPTAs("EMB", nontoxic_compts_EMB, ...
            AUCs_oral_EMB, AUCs_lung_EMB, ...
            Cmaxs_oral_EMB, Cmaxs_lung_EMB, ...
            MICs_TB_EMB, ...
            EMB_AUC_target, EMB_Cmax_target, ...
            oral_dose_EMB, oral_dose_freq_EMB, ...
            lung_dose_EMB, lung_dose_freq_EMB, ...
            n_pts_EMB, n_days_EMB);
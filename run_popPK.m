%% Run PopPK
% Script to run the population pharamcokinetic model for all drugs in the
% study


run_RIF = false;
run_PZA = true;
run_INH = true;

n_pts = 100;
n_days = 1;
relevant_compts = {"Lung", "Liver"};

addpath("Methods\")


%% Rifampin (RIF)

if run_RIF
    drug_RIF = "RIF";
    toxic_compts_RIF = ["Liver", "Kidney"];
    
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
    
    AUC_target_RIF  = 271; % AUC/MIC, from Jayaram et al.
    Cmax_target_RIF = 175; % Cmax/MIC, from Gumbo et al.
    
    oral_dose_RIF = 600; % mg
    lung_dose_RIF = 600; % mg
    oral_dose_freq_RIF = 1; % doses/day
    lung_dose_freq_RIF = 1; % doses/day
    
    setUpPopPK( ...
        drug_RIF, relevant_compts, toxic_compts_RIF, ...
        MICs_TB_RIF, AUC_target_RIF, Cmax_target_RIF, ...
        n_pts, n_days, ...
        oral_dose_RIF, lung_dose_RIF, oral_dose_freq_RIF, lung_dose_freq_RIF ...
        )
end


%% Pyrazinamide (PZA)

if run_PZA
    drug_PZA = "PZA";
    toxic_compts_PZA = ["Liver", "Kidney"];
    
    MICs_TB_PZA =  {[0.031, 28]; % TODO
                    [0.062, 47];
                    [0.120, 51];
                    [0.250, 76];
                    [0.500, 81];
                    [1.000, 68];
                    [2.000, 0];
                    [4.000, 5];
                    [8.000, 0];
                    [16.00, 0]};
    
    AUC_target_PZA  = 120; % AUC/MIC, from Gumbo et al.
    Cmax_target_PZA = 0;
    
    oral_dose_PZA = 2000; % mg
    lung_dose_PZA = 2000; % mg
    oral_dose_freq_PZA = 1; % doses/day
    lung_dose_freq_PZA = 1; % doses/day
    
    setUpPopPK( ...
        drug_PZA, relevant_compts, toxic_compts_PZA, ...
        MICs_TB_PZA, AUC_target_PZA, Cmax_target_PZA, ...
        n_pts, n_days, ...
        oral_dose_PZA, lung_dose_PZA, oral_dose_freq_PZA, lung_dose_freq_PZA ...
        )
end


%% Isoniazid (INH), Fast Acetylator (FA) and Slow Acetylator (SA)

if run_INH
    toxic_compts_INH = ["Liver", "Brain"];
    
    MICs_TB_INH =  {[0.031, 28]; % TODO
                    [0.062, 47];
                    [0.120, 51];
                    [0.250, 76];
                    [0.500, 81];
                    [1.000, 68];
                    [2.000, 0];
                    [4.000, 5];
                    [8.000, 0];
                    [16.00, 0]};
    
    AUC_target_INH  = 567; % AUC/MIC, Alffenaar et al.
    Cmax_target_INH = 0;
    
    oral_dose_INH = 300; % mg
    lung_dose_INH = 300; % mg
    oral_dose_freq_INH = 1; % doses/day
    lung_dose_freq_INH = 1; % doses/day
    
    % FA
    drug_INH = "INHFA";
    setUpPopPK( ...
        drug_INH, relevant_compts, toxic_compts_INH, ...
        MICs_TB_INH, AUC_target_INH, Cmax_target_INH, ...
        n_pts, n_days, ...
        oral_dose_INH, lung_dose_INH, oral_dose_freq_INH, lung_dose_freq_INH ...
        )
    
    % SA
    drug_INH = "INHSA";
    setUpPopPK( ...
        drug_INH, relevant_compts, toxic_compts_INH, ...
        MICs_TB_INH, AUC_target_INH, Cmax_target_INH, ...
        n_pts, n_days, ...
        oral_dose_INH, lung_dose_INH, oral_dose_freq_INH, lung_dose_freq_INH ...
        )
end
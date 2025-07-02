function [vol_PDs, flow_PDs, flow_frac_PDs] = getParamPDs(BW)
% function to set up probability distributions of physiological parameters
% for later sampling; flow rates are handled separately based on whether or
% not they are based on 

%% Volumes
v_cv = 0.2; % coeff. of variation, from Lyons et al.

v_v_mean = 3.6; % venous; ind. of BW
v_a_mean = 1.8; % arterial; ind. of BW
v_lu_mean = 0.0076 * BW;
v_brain_mean = 0.02 * BW;
v_heart_mean = 0.0047 * BW;
v_adipose_mean = (0.2142 * BW) / 0.916;
v_muscle_mean = 0.4 * BW;
v_skin_mean = 0.0371 * BW;
v_kidney_mean = 0.0044 * BW;
v_bone_mean = (0.1429 * BW) / 1.92;
v_spleen_mean = 0.0026 * BW;
v_gut_mean = 0.0171 * BW;
v_liver_mean = 0.0257 * BW;
v_others_mean = 0.04264 * BW;
v_ln_mean = 0.274; % ind. of BW
v_pl_mean = 0.3* BW / 1000; % mL/kg to L

v_means =  [v_v_mean, v_a_mean, v_lu_mean, v_brain_mean, v_heart_mean, ...
            v_adipose_mean, v_muscle_mean, v_skin_mean, v_kidney_mean, ...
            v_bone_mean, v_spleen_mean, v_gut_mean, v_liver_mean, ...
            v_others_mean, v_ln_mean, v_pl_mean];

v_sds = v_means .* v_cv;

% create cell array of truncated probability dists for each vol param
vol_PDs = cell(1, length(v_means));
for param_idx = 1:length(vol_PDs)
    current_mean = v_means(param_idx);
    current_sd = v_sds(param_idx);

    current_pd = makedist("Normal", "mu", current_mean, "sigma", current_sd);
    current_pd = truncate(current_pd, current_mean - 3 * current_sd, current_mean + 3 * current_sd);

    vol_PDs{param_idx} = current_pd;
end


%% Flows
q_cv = 0.3; % coeff. of variation, from Lyons et al.

q_pl_mean = 0.15 * BW/1000; % Pleural fluid flow [L/h]
q_bELF_mean = (5.25 * 100) * (5.75 * 10^-5) * (1/1000) * 3600; % from Himstedt et al.
q_aELF_mean = (171 * 100) * (5.75 * 10^-5) * (1/1000) * 3600; % from Himstedt et al.

q_la_mean_frac = 0.06;
q_sp_mean_frac = 77/5200; % spleen flow in L/h
q_gu_mean_frac = 1100/5200; % gut flow in L/h
q_br_mean_frac = 0.12;
q_hr_mean_frac = 0.04;
q_ad_mean_frac = 0.05;
q_mu_mean_frac = 0.17;
q_bo_mean_frac = 0.05;
q_sk_mean_frac = 0.05;
q_oth_mean_frac = 0.04365;
q_kd_mean_frac = 0.19;

% handle fracs and actual means separately for later processing
q_means = [q_pl_mean, q_bELF_mean, q_aELF_mean];

q_means_fracs = [q_la_mean_frac, q_sp_mean_frac, q_gu_mean_frac, ...
                q_br_mean_frac, q_hr_mean_frac, q_ad_mean_frac, ...
                q_mu_mean_frac, q_bo_mean_frac, q_sk_mean_frac, ...
                q_oth_mean_frac, q_kd_mean_frac];

q_sds = q_means .* q_cv;
q_sds_fracs = q_means_fracs .* q_cv;

% create cell arrays of truncated probability dists for each flow param
flow_PDs = cell(1, length(q_means));
for param_idx = 1:length(flow_PDs)
    current_mean = q_means(param_idx);
    current_sd = q_sds(param_idx);

    current_pd = makedist("Normal", "mu", current_mean, "sigma", current_sd);
    current_pd = truncate(current_pd, current_mean - 3 * current_sd, current_mean + 3 * current_sd);

    flow_PDs{param_idx} = current_pd;
end

flow_frac_PDs = cell(1, length(q_means_fracs));
for param_idx = 1:length(flow_frac_PDs)
    current_mean = q_means_fracs(param_idx);
    current_sd = q_sds_fracs(param_idx);

    current_pd = makedist("Normal", "mu", current_mean, "sigma", current_sd);
    current_pd = truncate(current_pd, current_mean - 3 * current_sd, current_mean + 3 * current_sd);

    flow_frac_PDs{param_idx} = current_pd;
end


end
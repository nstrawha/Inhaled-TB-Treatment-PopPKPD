% test script to analyze distributions of params based on AUC percentiles
% to be run after run_pop_PK_RIF.m with one relevant compt

%% TODO: Options for analysis

plot_hists = false;         % plot hists of AUCs and Cmaxs
params_regression = false;  % run regression for AUCs for all params
sum_effsize = true;         % run regression for the sum of all volume and flow eff sizes

r_cutoff = 0.3; % correlation coeff. past which to plot regressions

%% Helper function
% convert struct array to numeric matrix by extracting fields in fieldNames
function mat = struct2mat(structArray, fieldNames)
    n = numel(structArray);
    m = numel(fieldNames);
    mat = zeros(n, m);
    for i = 1:n
        for j = 1:m
            mat(i,j) = structArray(i).(fieldNames{j});
        end
    end
end


%% Plot hists of metrics

if plot_hists
    oral_AUCs = AUCs_oral_RIF{1, 1};
    lung_AUCs = AUCs_lung_RIF{1, 1};
    oral_Cmaxs = Cmaxs_oral_RIF{1, 1};
    lung_Cmaxs = Cmaxs_lung_RIF{1, 1};
    
    metrics = {oral_AUCs, lung_AUCs, oral_Cmaxs, lung_Cmaxs};
    titles = ["Oral Dose AUCs", "Lung Dose AUCs", "Oral Dose Cmaxs", "Lung Dose Cmaxs"];
    
    num_bins = floor(sqrt(length(oral_AUCs)));
    if num_bins > 30
        num_bins = 30;
    end
    
    for idx = 1:length(metrics)
        figure()
        histogram(metrics{idx}, num_bins)
        title(titles(idx))
    end

end


%% Unpack and cluster data

AUCs = AUCs_lung_RIF{1, 1};

% extract physiological params
phys_params = cellfun(@(x) x(8), param_store_RIF);

v_params = cellfun(@(x) x.V, phys_params);
q_params = cellfun(@(x) x.Q, phys_params);


%% Sort params into nested cell arrays for downstream analysis

% define field names for volumes and flows
raw_V_fields = {'V', 'A', 'LN'};
frac_V_fields = {'Lu', 'Brain', 'Heart', 'Adipose', 'Muscle', 'Skin', 'Kidney', ...
                 'Bone', 'Spleen', 'Gut', 'Liver', 'Others', 'Pl'};
raw_Q_fields = {'Pl', 'bELF', 'aELF'};
frac_Q_fields = {'LA', 'Sp', 'Gu', 'Br', 'Hr', 'Ad', 'Mu', ...
                 'Bo', 'Sk', 'Oth', 'Kd'};

% restructure volume params storage
raw_structs = arrayfun(@(s) rmfield(s, setdiff(fieldnames(s), raw_V_fields)), v_params);
frac_structs = arrayfun(@(s) rmfield(s, setdiff(fieldnames(s), frac_V_fields)), v_params);
bw_vec = (arrayfun(@(s) s.BW, v_params))';
v_fsum_vec = (arrayfun(@(s) s.fsum, v_params))';

% convert struct arrays to numeric matrices
raw_v_matrix = struct2mat(raw_structs, raw_V_fields);
frac_v_matrix = struct2mat(frac_structs, frac_V_fields);

% backcalculate fractional fields by BW and fsum
frac_v_matrix = (frac_v_matrix ./ bw_vec) .* v_fsum_vec;

% store results
v_params = {bw_vec, raw_v_matrix, frac_v_matrix};
v_params_mat = [bw_vec, raw_v_matrix, frac_v_matrix];

% restructure flow params storage
raw_structs = arrayfun(@(s) rmfield(s, setdiff(fieldnames(s), raw_Q_fields)), q_params);
frac_structs = arrayfun(@(s) rmfield(s, setdiff(fieldnames(s), frac_Q_fields)), q_params);
qc_vec = (arrayfun(@(s) s.QC, q_params))';
q_fsum_vec = (arrayfun(@(s) s.fsum, q_params))';

% convert struct arrays to numeric matrices
raw_q_matrix = struct2mat(raw_structs, raw_Q_fields);
frac_q_matrix = struct2mat(frac_structs, frac_Q_fields);

% backcalculate fractional fields by BW and fsum
frac_q_matrix = (frac_q_matrix ./ qc_vec) .* q_fsum_vec;

% store results
q_params = {qc_vec, raw_q_matrix, frac_q_matrix};
q_params_mat = [qc_vec, raw_q_matrix, frac_q_matrix];

% store all params in nested cell array
all_params = {v_params, q_params};
all_params_mats = {v_params_mat, q_params_mat};


%% Package PDs

% package PDs into nested cell array corresponding to param storage
% UNUSED
v_PDs = {{bw_PD}, vol_PDs_RIF, vol_frac_PDs_RIF};
q_PDs = {{qc_PD}, flow_PDs_RIF, flow_frac_PDs_RIF};
v_PDs = [v_PDs{:}];
q_PDs = [q_PDs{:}];

all_PDs = {v_PDs, q_PDs};

% set up names
type_names = {"Volume", "Flow"};
class_names = {"BW/QC", "Raw", "Fractional"};
param_names = {{{"Body Weight"}, ... % volumes
                {"Venous Blood", "Arterial Blood", "Lymph Node"}, ...
                {"Lung", "Brain", "Heart", "Adipose", "Muscle", "Skin", "Kidney", "Bone", "Spleen", "Gut", "Liver", "Others", "Pleura"}}, ...
               {{"Cardiac Output"}, ... % flows
                {"Pleura", "bELF", "aELF"}, ...
                {"LA", "Spleen", "Gut", "Brain", "Heart", "Adipose", "Muscle", "Bone", "Skin", "Other", "Kidney", "Liver"}}};


%% Conduct analysis
% iterate through param types (vols and flows)
for type_idx = 1:length(all_params)
    current_type = all_params{type_idx};

    type_name = type_names{type_idx};
    current_type_names = param_names{type_idx};

    current_params_mat = all_params_mats{type_idx};

    % iterate through param classes (BW/QC, raw, frac)
    for class_idx = 1:length(current_type)
        current_param_set = current_type{class_idx};

        class_name = class_names{class_idx};
        current_class_names = current_type_names{class_idx};

        % iterate through params 
        for param_idx = 1:size(current_param_set, 2) % ncol = nparams
            current_params = current_param_set(:, param_idx);

            param_name = current_class_names{param_idx};

            % carry out regression with all params
            if params_regression

                % make scatterplot of param values vs AUC if r > 0.5
                corr_coeff = corrcoef(current_params, lung_AUCs);
    
                if abs(corr_coeff(1,2)) >= r_cutoff
                    figure()
                    scatter(current_params, lung_Cmaxs, ...
                        "MarkerEdgeColor", "none", ...
                        "MarkerFaceColor", "black")
                    lsline;
   
                    xlabel("Parameter Value")
                    ylabel("AUC")
                    title(sprintf('Param Type: %s; Class: %s; Param: %s; r = %.3f', ...
                                    type_name, class_name, param_name, corr_coeff(1,2)));
                end
            end
            
        end
    end

    % carry out regression with sum of effect sizes
    if sum_effsize
        current_type_pds = all_PDs{type_idx};
        sum_effsizes_store = zeros(1, length(lung_AUCs));

        % iterate through rows (patients)
        for pt_idx = 1:size(current_params_mat, 1)
            current_row = current_params_mat(pt_idx, :);

            % iterate through individual parameters
            for param_idx = 1:length(current_row)
                current_param = current_row(param_idx);
                current_pd = current_type_pds{param_idx};

                pdmean = mean(current_pd);
                pdsd = std(current_pd);

                % calculate and store effect size
                effsize = (current_param - pdmean) / pdsd;
                sum_effsizes_store(pt_idx) = sum_effsizes_store(pt_idx) + effsize;
            end
        end

        % create scatterplot of results
        corr_coeff = corrcoef(sum_effsizes_store, lung_AUCs);

        figure()
        scatter(sum_effsizes_store, lung_AUCs, ...
            "MarkerEdgeColor", "none", ...
            "MarkerFaceColor", "black")
        lsline;

        xlabel("Sum of Parameter Effect Sizes")
        ylabel("AUC")
        title(sprintf('Param Type: %s; r = %.3f', ...
                        type_name, corr_coeff(1,2)))
    end

end

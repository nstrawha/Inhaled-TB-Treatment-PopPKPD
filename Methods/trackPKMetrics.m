function trackPKMetrics(compt_list, params, toxic_compts, odose, odose_freq, ldose, ldose_freq, compts_to_plot, days, t, Cs_oral, Cs_lung)
% function to analyze AUC_24, C_avg, and C_max and 
% write results to output file

%% Set up cell array formatting
cells_store = cell(1, length(compts_to_plot));
tbl_format = cell(11, 4);

tbl_format{1, 1} = append("Day ", num2str(days));

tbl_format{2, 1} = "Oral Dose";
tbl_format{3, 1} = append(num2str(odose), " mg, ", num2str(odose_freq), " x/day");
tbl_format{5, 1} = "Lung Dose";
tbl_format{6, 1} = append(num2str(ldose), " mg, ", num2str(ldose_freq), " x/day");
tbl_format{9, 1} = "Better Dose";

tbl_format{1, 2} = "Metric Type";
tbl_format{2, 2} = "AUC_24";
tbl_format{3, 2} = "C_avg (ug/mL)";
tbl_format{4, 2} = "C_max (ug/mL)";
tbl_format{5, 2} = "AUC_24";
tbl_format{6, 2} = "C_avg, (ug/mL)";
tbl_format{7, 2} = "C_max (ug/mL)";
tbl_format{9, 2} = "AUC_24";
tbl_format{10, 2} = "C_avg (ug/mL)";
tbl_format{11, 2} = "C_max (ug/mL)";

tbl_format{1, 3} = "Mean";
tbl_format{1, 4} = "SD";
tbl_format{8, 3} = "p-value";
tbl_format{8, 4} = "Effect Size";

last_day_start = length(t) - 24 / params{14} + 1;


%% Iterate through compartments
for compt_idx = 1:length(compts_to_plot)

    current_compt = compts_to_plot{compt_idx};
    current_cell = tbl_format;
    [idx_to_calc, ~] = find(string(compt_list) == current_compt);

    % pull compt specific timecourses
    current_cs_oral = cell2mat(cellfun(@(x) x(:, idx_to_calc), Cs_oral, "UniformOutput", false));
    current_cs_lung = cell2mat(cellfun(@(x) x(:, idx_to_calc), Cs_lung, "UniformOutput", false));

    % calculate metrics
    AUCs_oral = trapz(t(last_day_start:end), current_cs_oral(last_day_start:end, :));
    AUCs_lung = trapz(t(last_day_start:end), current_cs_lung(last_day_start:end, :));

    Cavgs_oral = mean(current_cs_oral(last_day_start:end, :));
    Cavgs_lung = mean(current_cs_lung(last_day_start:end, :));

    Cmaxs_oral = max(current_cs_oral(last_day_start:end, :));
    Cmaxs_lung = max(current_cs_lung(last_day_start:end, :));

    % calculate metric info for comparison
    AUCs_oral_mean  = mean(AUCs_oral);
    AUCs_oral_sd    = std(AUCs_oral);
    AUCs_lung_mean  = mean(AUCs_lung);
    AUCs_lung_sd    = std(AUCs_lung);

    Cavgs_oral_mean  = mean(Cavgs_oral);
    Cavgs_oral_sd    = std(Cavgs_oral);
    Cavgs_lung_mean  = mean(Cavgs_lung);
    Cavgs_lung_sd    = std(Cavgs_lung);

    Cmaxs_oral_mean  = mean(Cmaxs_oral);
    Cmaxs_oral_sd    = std(Cmaxs_oral);
    Cmaxs_lung_mean  = mean(Cmaxs_lung);
    Cmaxs_lung_sd    = std(Cmaxs_lung);

    % calculate better method for each metric
    if ~ismember(current_compt, toxic_compts)

        if AUCs_oral_mean >= AUCs_lung_mean
            better_AUC = "Oral";
        else
            better_AUC = "Lung";
        end

        if Cavgs_oral_mean >= Cavgs_lung_mean
            better_Cavg = "Oral";
        else
            better_Cavg = "Lung";
        end

        if Cmaxs_oral_mean >= Cmaxs_lung_mean
            better_Cmax = "Oral";
        else
            better_Cmax = "Lung";
        end

    else
        if AUCs_oral_mean > AUCs_lung_mean
            better_AUC = "Lung";
        else
            better_AUC = "Oral";
        end

        if Cavgs_oral_mean > Cavgs_lung_mean
            better_Cavg = "Lung";
        else
            better_Cavg = "Oral";
        end

        if Cmaxs_oral_mean > Cmaxs_lung_mean
            better_Cmax = "Lung";
        else
            better_Cmax = "Oral";
        end
    end

    % calculate effect sizes 
    AUCs_effsize  = meanEffectSize(AUCs_oral, AUCs_lung);
    Cavgs_effsize = meanEffectSize(Cavgs_oral, Cavgs_lung);
    Cmaxs_effsize = meanEffectSize(Cmaxs_oral, Cmaxs_lung);

    % perform hypothesis tests
    [~, AUC_p]  = ttest(AUCs_oral, AUCs_lung);
    [~, Cavg_p] = ttest(Cavgs_oral, Cavgs_lung);
    [~, Cmax_p] = ttest(Cmaxs_oral, Cmaxs_lung);

    % record info
    % means
    current_cell{2, 3} = round(AUCs_oral_mean, 2);
    current_cell{3, 3} = round(Cavgs_oral_mean, 2);
    current_cell{4, 3} = round(Cmaxs_oral_mean, 2);
    current_cell{5, 3} = round(AUCs_lung_mean, 2);
    current_cell{6, 3} = round(Cavgs_lung_mean, 2);
    current_cell{7, 3} = round(Cmaxs_lung_mean, 2);

    % SDs
    current_cell{2, 4} = round(AUCs_oral_sd, 2);
    current_cell{3, 4} = round(Cavgs_oral_sd, 2);
    current_cell{4, 4} = round(Cmaxs_oral_sd, 2);
    current_cell{5, 4} = round(AUCs_lung_sd, 2);
    current_cell{6, 4} = round(Cavgs_lung_sd, 2);
    current_cell{7, 4} = round(Cmaxs_lung_sd, 2);

    % tests
    current_cell{9, 3}  = append(better_AUC, ", p = ", num2str(AUC_p));
    current_cell{10, 3} = append(better_Cavg, ", p = ", num2str(Cavg_p));
    current_cell{11, 3} = append(better_Cmax, ", p = ", num2str(Cmax_p, 3));

    current_cell{9, 4}  = abs(round(AUCs_effsize.Effect, 2));
    current_cell{10, 4} = abs(round(Cavgs_effsize.Effect, 2));
    current_cell{11, 4} = abs(round(Cmaxs_effsize.Effect, 2));

    cells_store{compt_idx} = current_cell;

end


%% Write output file
for compt_idx = 1:length(compts_to_plot)
    writecell(cells_store{compt_idx}, append("Outputs/popPK_analysis_day", num2str(days), ".xlsx"), "Sheet", compts_to_plot{compt_idx});
end


end
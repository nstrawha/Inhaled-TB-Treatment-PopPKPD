function plotPTAs(drug, tval_type, nontoxic_compts, AUCs_oral, AUCs_lung, Cmaxs_oral, Cmaxs_lung, MIC_dist, target, oral_dose, oral_dose_freq, lung_dose, lung_dose_freq, n_pts, n_days)
% PLOTPTAS - Calculates and plots the probability of target attainment (for
% AUC/MIC and Cmax/MIC targets) for an oral vs. lung dose of medication.
%
% DESCRIPTION:
%   This function uses the PK metrics calculated for all patients and a
%   drug-TB MIC distribution in order to calculate and plot the probability 
%   of target attainment (PTA) for each metric in each compartment.
%
% INPUTS:
% - drug (str): An all-caps three-letter identifier of the relevant drug
% - tval_type (str): An identifier for the type of target value to
%   calculate
% - nontoxic compts (str array): Contains the names of all nontoxic
%   compartments for which PTA is to be calculated and plotted
% - AUCs_oral (cell array): Contains AUCs for all patients given an 
%   oral dose in each cell, with each cell representing a different 
%   compartment
% - AUCs_lung (cell array): Contains AUCs for all patients given a 
%   lung dose in each cell, with each cell representing a different 
%   compartment
% - Cmaxs_oral (cell array): Contains Cmaxs for all patients given an 
%   oral dose in each cell, with each cell representing a different 
%   compartment
% - Cmaxs_lung (cell array): Contains Cmaxs for all patients given a
%   lung dose in each cell, with each cell representing a different 
%   compartment
% - MIC_dist (cell array): Contains arrays formatted like [MIC, # isolates]
%   describing the drug-TB MIC distribution
% - target (double): The AUC/MIC target for the drug
% - target (double): The Cmax/MIC target for the drug
% - oral_dose (int): The oral dose amt. of the drug (in mg)
% - oral_dose_freq (int): The number of times per day an oral dose is to be
%   administered
% - lung_dose (int): The lung dose amt. of the drug (in mg)
% - lung_dose_freq (int): The number of times per day an inhaled/lung dose
%   is to be administered
% - n_pts (int): The number of patients for which concentration-time
%   courses have been calculated
% - n_days (int): The number of days for which concentration-time courses
%   have been calculated

% unpack MIC dist info
concs = cell2mat(cellfun(@(x) x(:, 1), MIC_dist, "UniformOutput", false));
isolates = cell2mat(cellfun(@(x) x(:, 2), MIC_dist, "UniformOutput", false));

MICs_to_test = 0.001:0.001:25; % start at 0.001 ug/mL

% set up storage
PTA_storage = cell(1, length(nontoxic_compts));
CFR_storage = cell(1, length(nontoxic_compts));

% find the continuous lognormal probabiliity dist for MICs
dummy_data = repelem(concs, isolates);
MIC_PD = fitdist(dummy_data, "Lognormal");


%% Calculate PTA dists

% AUC/MIC tval type
if tval_type == "AUC/MIC"
    parfor compt_idx = 1:length(nontoxic_compts)
    
        current_AUCs_oral = AUCs_oral{compt_idx};
        current_AUCs_lung = AUCs_lung{compt_idx};
    
        PTA_array = zeros(2, length(MICs_to_test));
    
        % iterate through MICs
        for mic_idx = 1:length(MICs_to_test)
    
            current_mic = MICs_to_test(mic_idx);
    
            success_pts_AUC_oral = 0;
            success_pts_AUC_lung = 0;
    
            % iterate through patients
            for pt_idx = 1:n_pts
                pt_AUC_oral = current_AUCs_oral(pt_idx);
                pt_AUC_lung = current_AUCs_lung(pt_idx);
    
                % track AUC/MIC successes
                if pt_AUC_oral/current_mic >= target
                    success_pts_AUC_oral = success_pts_AUC_oral + 1;
                end
    
                if pt_AUC_lung/current_mic >= target
                    success_pts_AUC_lung = success_pts_AUC_lung + 1;
                end
    
            end
    
            PTA_array(1:2, mic_idx) = [success_pts_AUC_oral / n_pts; success_pts_AUC_lung / n_pts];
    
        end
    
        PTA_storage{compt_idx} = PTA_array;
    
    end

% Cmax/MIC tval type
elseif tval_type == "Cmax/MIC"
        parfor compt_idx = 1:length(nontoxic_compts)
    
        current_Cmaxs_oral = Cmaxs_oral{compt_idx};
        current_Cmaxs_lung = Cmaxs_lung{compt_idx};

        PTA_array = zeros(2, length(MICs_to_test));
    
        % iterate through MICs
        for mic_idx = 1:length(MICs_to_test)
    
            current_mic = MICs_to_test(mic_idx);
    
            success_pts_Cmax_oral = 0;
            success_pts_Cmax_lung = 0;
    
            % iterate through patients
            for pt_idx = 1:n_pts
                pt_Cmax_oral = current_Cmaxs_oral(pt_idx);
                pt_Cmax_lung = current_Cmaxs_lung(pt_idx);
    
                % track Cmax/MIC successes
                if pt_Cmax_oral/current_mic >= target
                    success_pts_Cmax_oral = success_pts_Cmax_oral + 1;
                end
    
                if pt_Cmax_lung/current_mic >= target
                    success_pts_Cmax_lung = success_pts_Cmax_lung + 1;
                end
    
            end
    
            PTA_array(1:2, mic_idx) = [success_pts_Cmax_oral / n_pts; success_pts_Cmax_lung / n_pts];
    
        end
    
        PTA_storage{compt_idx} = PTA_array;
    
    end

end

%% Plot PTA results

for compt_idx = 1:length(PTA_storage)

    current_PTA_array = PTA_storage{compt_idx};
    oral_PTAs = current_PTA_array(1, :);
    lung_PTAs = current_PTA_array(2, :);

    fig = figure();
    hold on;

    % MIC distribution
    plot(MICs_to_test, pdf(MIC_PD, MICs_to_test) / max(pdf(MIC_PD, MICs_to_test)), ...
        "DisplayName", "Continuous MIC Dist. Estimate (Normalized)", ...
        "LineWidth", 2, ...
        "Color", "Black");

    % PTAs
    plot(MICs_to_test, oral_PTAs, ...
        "DisplayName", append("Oral PTA; ", num2str(oral_dose), " mg, ", num2str(oral_dose_freq), " x/day"), ...
        "LineWidth", 1.5, "Color", "Blue");
    plot(MICs_to_test, lung_PTAs, ...
        "DisplayName", append("Lung PTA; ", num2str(lung_dose), " mg, ", num2str(lung_dose_freq), " x/day"), ...
        "LineWidth", 1.5, "Color", "Red");
    
    yline(0, "LineWidth", 2, "HandleVisibility", "off");
    xlim([min(MICs_to_test), max(MICs_to_test)]);
    ylim([-0.05, 1.05]);
    xlabel("Minimum Inhibitory Concentration (\mug/mL)");
    ylabel("Probability of Target Attainment (PTA)");
    title(append(nontoxic_compts{compt_idx}, " PTA ", tval_type), "FontSize", 20);

    set(gca, "XScale", "log");
    leg = legend("show");
    leg.FontSize = 8;
    hold off;
    
    % predict CFR
    oral_inhibition_yvals = oral_PTAs .* pdf(MIC_PD, MICs_to_test);
    lung_inhibition_yvals = lung_PTAs .* pdf(MIC_PD, MICs_to_test);

    CFR_oral = trapz(MICs_to_test, oral_inhibition_yvals);
    CFR_lung = trapz(MICs_to_test, lung_inhibition_yvals);

    CFR_storage{compt_idx} = [CFR_oral, CFR_lung];

    % save figure
    saveas(fig, append("Outputs/", drug, "/Figures/", nontoxic_compts{compt_idx}, "_", erase(tval_type, "/"), ".png"));

end


%% CFR barplots

for compt_idx = 1:length(PTA_storage)
    current_compt = nontoxic_compts{compt_idx};
    
    fig = figure();

    current_bar = bar(1:2, [CFR_oral * 100, CFR_lung * 100], 0.5);
    current_bar.FaceColor = "flat";
    current_bar.CData(1, :) = [0 0 1];
    current_bar.CData(2, :) = [1 0 0];

    bar_heights = [CFR_oral * 100, CFR_lung * 100];
    label_text1 = append(num2str(round(CFR_oral * 100, 2)), "%");
    label_text2 = append(num2str(round(CFR_lung * 100, 2)), "%");
    text(1, bar_heights(1) + 2, label_text1, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');
    text(2, bar_heights(2) + 2, label_text2, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');

    ylim([0, 100]);
    xticks(1:2)
    xticklabels(["Oral Dose", "Lung Dose"]);
    title(append("CFR for ", tval_type, " Target"), "FontSize", 15);
    ylabel("% CFR");

    saveas(fig, append("Outputs/", drug, "/Figures/", current_compt, "_CFRs_comparison.png"));

end


end
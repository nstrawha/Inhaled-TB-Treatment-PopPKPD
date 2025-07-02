function plotTimecourses(params, compts_to_plot, ndays, days_to_plot, oral_dose_freq, lung_dose_freq, t_oraldose, Cs_oraldose, t_lungdose, Cs_lungdose)
% function to plot timecourses in a certain compartment (or multiple
% compartments) for comparison between oral and lung dosing

% check for input error
if days_to_plot > ndays
    disp("Error: cannot plot more days than days calculated for")
    return
end

% set up figure as: oral conc | lung conc
figure();
tiledlayout(1, 2);
set(0, "DefaultFigureWindowStyle", "docked");

%% Plotting
% oral dose
nexttile
plot(t_oraldose, Cs_oraldose, "LineWidth", 0.5);

xlabel('Time (h)'); ylabel('Concentration (\mug/mL)');
xlim([(ndays - 1 * days_to_plot) * 24, ndays * 24]);
ylim([0, 40]);
title(append("Oral Dose (", num2str(params{1}), " mg, ", num2str(oral_dose_freq), "x/day)"));
set(gca,'FontSize', 20);
grid on;

% lung dose
nexttile
plot(t_oraldose, Cs_lungdose, "LineWidth", 0.5);

xlabel('Time (h)'); ylabel('Concentration (\mug/mL)');
xlim([(ndays - 1 * days_to_plot) * 24, ndays * 24]);
ylim([0, 40]);
title(append("Lung Dose (", num2str(params{2}), " mg, ", num2str(lung_dose_freq), "x/day)"));
set(gca,'FontSize', 20);
grid on;

if ndays <= 7
    if isscalar(compts_to_plot)
        sgtitle(append(compts_to_plot, " Concentration, Day ", num2str(ndays)), "FontSize", 30);
    
    else
        sgtitle(append("Concentrations in Various Sites, Day ", num2str(ndays)), "FontSize", 25);
        legend show;
    
    end
else
    if isscalar(compts_to_plot)
        sgtitle(append(compts_to_plot, " Concentration, Steady-State "), "FontSize", 30);
    
    else
        sgtitle(append("Concentrations in Various Sites, Steady-State "), "FontSize", 25);
        legend show;
    
    end
end


end
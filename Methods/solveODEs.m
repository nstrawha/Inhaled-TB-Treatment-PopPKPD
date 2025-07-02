function [t_oraldose, C_oraldose, t_lungdose, C_lungdose] = solveODEs(drug, params, ncompts_total, ndays, oral_dose_freq, lung_dose_freq)
% function to solve ODEs for oral vs. lung dosing; ODEs from
% Ramachandran & Gadgil, 2023

options = odeset('RelTol',1e-6,'AbsTol',1e-8);

%% Unpack parameters
if drug == "RIF"
    oral_dose = params{1};
    lung_dose = params{2};
    ka_oral = params{3};
    CL = params{4};
    fR = params{5};
    kr = params{6};
    pt = params{7};
    phys = params{8};
    kF = params{9};
    ka_lung = params{10};
    effRB = params{11};
    effRA = params{12};
    br_frac = params{13};
    tstep = params{14};

else
    disp("Invalid drug specified");
    return

end

%% Set up matrices to track solutions
C0_oraldose = zeros(1, ncompts_total); C0_oraldose(end) = oral_dose;
C0_lungdose = zeros(1, ncompts_total); C0_lungdose(end) = lung_dose;

timepts_oral = 0:tstep:(24 / oral_dose_freq);
timepts_lung = 0:tstep:(24 / lung_dose_freq);
t_oraldose = 0:tstep:(24 * ndays - tstep);
t_lungdose = 0:tstep:(24 * ndays - tstep);

C_oraldose = zeros((length(timepts_oral) - 1) * oral_dose_freq * ndays, ncompts_total);
C_lungdose = zeros((length(timepts_lung) - 1) * lung_dose_freq * ndays, ncompts_total);

%% Solve ODEs
if drug == "RIF"
    % oral eqs
    for dose_idx = 1:(ndays * oral_dose_freq)
        
        [~, C_oraldose_temp] = ode23s(@(t, A) RIF_oral_ODEs(t, A, ka_oral, kr, kF, ...
                            CL, fR, phys, pt), timepts_oral, C0_oraldose, options);
        C0_oraldose = [C_oraldose_temp(end, 1:(ncompts_total - 1))'; 
                                    C_oraldose_temp(end, ncompts_total) + oral_dose];
        C_oraldose_temp(end, :) = []; % remove initial condition

        % store result
        temp_row_start = (dose_idx - 1) * (length(timepts_oral) - 1) + 1;
        temp_row_end = dose_idx * (length(timepts_oral) - 1);
        C_oraldose(temp_row_start:temp_row_end, :) = C_oraldose_temp;
    end

    % lung eqs
    for dose_idx = 1:(ndays * lung_dose_freq)
        
        [~, C_lungdose_temp] = ode23s(@(t, A) RIF_lung_ODEs(t, A, ka_lung, kr, kF, effRB, ...
                            effRA, br_frac, CL, fR, phys, pt), timepts_lung, C0_lungdose, options);
        C0_lungdose = [C_lungdose_temp(end, 1:(ncompts_total - 1))'; 
                                    C_lungdose_temp(end, ncompts_total) + lung_dose];
        C_lungdose_temp(end, :) = []; % remove initial condition

        % store result
        temp_row_start = (dose_idx - 1) * (length(timepts_lung) - 1) + 1;
        temp_row_end = dose_idx * (length(timepts_lung) - 1);
        C_lungdose(temp_row_start:temp_row_end, :) = C_lungdose_temp;
    end

else
    disp("Invalid drug specified");
    return

end


end
function phys = loadPhysiology(BW)
% Load physiological values from Table S2 and S3
% BW = Body weight (kg)

% Volumes in Liters
phys.V.V = 3.6; % Venous
phys.V.A = 1.8; % Arterial
phys.V.Lu = 0.0076 * BW;
phys.V.Brain = 0.02 * BW;
phys.V.Heart = 0.0047 * BW;
phys.V.Adipose = 0.2142 * BW;
phys.V.Muscle = 0.4 * BW;
phys.V.Skin = 0.0371 * BW;
phys.V.Kidney = 0.0044 * BW;
phys.V.Bone = 0.1429 * BW;
phys.V.Spleen = 0.0026 * BW;
phys.V.Gut = 0.0171 * BW;
phys.V.Liver = 0.0257 * BW;
phys.V.Others = 0.04264 * BW;
phys.V.LN = 0.274; % Fixed volume
phys.V.Pl = 0.3 * BW / 1000; % mL/kg to L

% Cardiac output in L/h
QC = 5200 / 1000 * 60; % mL/min to L/h
phys.Q.total = QC;
phys.Q.LA = 0.06 * QC;
phys.Q.Sp = 77 / 1000 * 60; % spleen flow in L/h
phys.Q.Pl = 0.15 * BW/1000;  % Pleural fluid flow [L/h]
phys.Q.Gu = 1100 / 1000 * 60; % gut flow in L/h
% brain, adipose, heart,  muscle, skin, others, bone, spleen,kidney, gut,
% liver
phys.Q.Tissue = QC * [0.12,0.05,0.04,  0.17,  0.05,0.04365,...
    0.05,77/5200, ...
    0.19, 1100/5200, 0.06 + 1100/5200 + 77/5200];

% Lymph flow in L/h
lymph_total = 8 / 24; % L/day to L/h
phys.L.Lu = 0.03 * lymph_total;
% brain, adipose, heart,  muscle, skin, others, bone, spleen,kidney, gut,
% liver
phys.L.Tissue = lymph_total * [0.0105, 0.128,0.01,  0.16, 0.0703,0.0562,...
    0,0,...
    0.085, 0.12, 0.33];
phys.L.LN = lymph_total;
phys.L.Li = 0.33 * lymph_total;
end

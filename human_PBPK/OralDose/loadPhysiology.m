function phys = loadPhysiology(BW)
% Load physiological values from Table S2 and S3
% BW = Body weight (kg)

% Volumes in Liters
phys.V.V = 3.6; % Venous
phys.V.A = 1.8; % Arterial
phys.V.Lu = 0.0076 * BW;
phys.V.Brain = 0.02 * BW;
phys.V.Heart = 0.0047 * BW;
phys.V.Adipose = (0.2142 * BW)/0.916;
phys.V.Muscle = 0.4 * BW;
phys.V.Skin = 0.0371 * BW;
phys.V.Kidney = 0.0044 * BW;
phys.V.Bone = (0.1429 * BW)/1.92;
phys.V.Spleen = 0.0026 * BW;
phys.V.Gut = 0.0171 * BW;
phys.V.Liver = 0.0257 * BW;
phys.V.Others = 0.04264 * BW;
phys.V.LN = 0.274; % Fixed volume
phys.V.Pl = 0.3 * BW / 1000; % mL/kg to L

phys.V.Tissues = [phys.V.Brain, phys.V.Adipose, phys.V.Heart, phys.V.Muscle, phys.V.Skin, phys.V.Others];
% Cardiac output in L/h
QC = 5200 / 1000 * 60; % mL/min to L/h
phys.Q.total = QC;
phys.Q.LA = 0.06 * QC;
phys.Q.Sp = 77/1000 * 60 ; % spleen flow in L/h
phys.Q.Pl = 0.15 * BW/1000;  % Pleural fluid flow [L/h]
phys.Q.Gu = 1100/1000 *60; % gut flow in L/h

phys.Q.Br = QC * 0.12;
phys.Q.Hr = QC *0.04;
phys.Q.Ad = QC *0.05;
phys.Q.Mu = QC *0.17;
phys.Q.Bo = QC *0.05;
phys.Q.Sk = QC *0.05;
phys.Q.Oth = QC *0.04365;
phys.Q.Kd = QC *0.19;
phys.Q.Li = phys.Q.Sp + phys.Q.Gu + phys.Q.LA;

phys.Q.TissuesV = [phys.Q.Br, phys.Q.Ad, phys.Q.Hr, phys.Q.Mu, phys.Q.Sk, phys.Q.Oth, phys.Q.Bo, phys.Q.Kd, phys.Q.Li];

% phys.Q.TissuesA = [phys.Q.Br, phys.Q.Ad, phys.Q.Hr, phys.Q.Mu, phys.Q.Sk, phys.Q.Oth, phys.Q.Bo,  phys.Q.Sp, phys.Q.Kd, phys.Q.Gu, phys.Q.LA];


% Lymph flow in L/h
lymph_total = 8 / 24; % L/day to L/h
phys.L.Lu = lymph_total*0.03;

phys.L.Br = lymph_total * 0.0105;
phys.L.Hr = lymph_total *0.01;
phys.L.Ad = lymph_total *0.128;
phys.L.Mu = lymph_total *0.16;
phys.L.Bo = 0;
phys.L.Sk = lymph_total *0.0703;
phys.L.Oth = lymph_total *0.0562;
phys.L.Kd = lymph_total *0.085;
phys.L.Li = lymph_total *0.33;
phys.L.Sp = 0;
phys.L.Gu = lymph_total *0.12;

phys.L.TissuesV = [phys.L.Br, phys.L.Ad, phys.L.Hr,phys.L.Mu, phys.L.Sk, phys.L.Oth,phys.L.Bo, phys.L.Kd, phys.L.Li];
phys.L.TissuesL = [phys.L.Br, phys. L.Ad, phys.L.Hr,phys.L.Mu,phys.L.Sk, phys.L.Oth, phys.L.Kd, phys.L.Gu,phys.L.Li,];

phys.L.LN = lymph_total;

end

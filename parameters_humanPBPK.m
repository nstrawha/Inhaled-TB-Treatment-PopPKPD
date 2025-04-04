% Parameters.m
function params = parameters_humanPBPK()
% Cardiac output = 5200 (mL/min) convert to L/hr
params.Qc = (5200*60)/1000;

% Body weight (kg)
params.bodyweight = 70; 

% Blood flow rates mL/hr
params.Qbr = 0.12*params.Qc;% brain
params.Qad = 0.05*params.Qc;%adipose
params.Qhr = 0.04*params.Qc;%heart
params.Qmu = 0.17*params.Qc; % muscle
params.Qbo = 0.05*params.Qc;% bone
params.Qsk = 0.05*params.Qc;% skin
params.Qkd = 0.19*params.Qc;% kidney
params.Qsp = 77;% spleen (77/5200)*5200
params.Qgu = 1100;% gut (1100/5200)*5200
params.Qla = 0.06*params.Qc; %hepatic artery
params.Qli = params.Qla + params.Qgu + params.Qsp;% liver
params.Qoth = 0.04365*params.Qc;% others
params.Qpl = (0.15*params.bodyweight)/1000;% pleura mL/hr/kg convert to L/hr by multiplying by bodyweight and dividing by 1000 mL
%Sum of all compartments minus lung
params.Qt = params.Qbr + params.Qad + params.Qhr + params.Qmu + params.Qbo + params.Qsk + params.Qkd + params.Qsp + params.Qgu + params.Qli + params.Qla + params.Qoth + params.Qpl;

% Lymph flow rates (L/day)
params.alf = 8/24; % Afferent lymph flow rate = 8 L/day convert to L/hr
params.Llu = 0.03 * params.alf; % Lungs
params.Lbr = 0.0105 * params.alf; % Brain
params.Llu = 0.03*params.alf;% lungs
params.Lbr = 0.0105*params.alf;% brain
params.Lad = 0.128*params.alf;%adipose
params.Lhr = 0.01*params.alf;%heart
params.Lmu = 0.16*params.alf; % muscle
params.Lbo = 0;% bone
params.Lsk = 0.0703*params.alf;% skin
params.Lkd = 0.085*params.alf;% kidney
params.Lsp = 0;% spleen
params.Lgu = 0.12*params.alf;% gut
params.Lli = 0.33*params.alf;% liver
params.Loth = 0.0562*params.alf;% others
params.Lln =params.alf;

%Sum of all compartments minus lung
params.Lt =  params.Lbr + params.Lad + params.Lhr + params.Lmu + params.Lbo + params.Lsk + params.Lkd + params.Lsp + params.Lgu + params.Lli + params.Loth;

%% Volumes (L)
% As fraction of body weight, assume 70 kg man
% Assume mass = volume and density = 1 g/cm3
params.Vlu = 0.0076 * params.bodyweight; % Lung
params.Vbr = 0.02 * params.bodyweight; % Brain
params.Vad = (0.2142*params.bodyweight)/0.916;%adipose density = 0.916 g/cm3
params.Vhr = 0.0047*params.bodyweight;%heart
params.Vmu = 0.4*params.bodyweight; % muscle
params.Vbo = (0.1429*params.bodyweight)/1.92;% bone density = 1.92 g/cm3
params.Vsk = 0.0371*params.bodyweight;% skin
params.Vkd = 0.0044*params.bodyweight;% kidney
params.Vsp = 0.0026*params.bodyweight;% spleen
params.Vgu = 0.0171*params.bodyweight;% gut
params.Vli = 0.0257*params.bodyweight;% liver
params.Vln = 0.274;% lymph node (0.274/70) * 70
params.Va = 1.8;% arterial blood (1.8/70) * 70
params.Vv = 3.6;% venous blood (3.6/70) * 70
params.Voth = 0.04264*params.bodyweight;% others
params.Vpl = (0.3*params.bodyweight)/1000;% pleura

%% Partition coefficients - calcuated with Rowland Rodgers equations
% From Ramachandran 2023 who calculated using RR equations Table S5
%For RIF

params.Plu = 1.7115; % Lung
params.Pbr = 0.2285; % Brain
params.Pad = 0.1885;%adipose
params.Phr = 1.0158;%heart
params.Pmu = 0.6949; % muscle
params.Pbo = 0.3157;% bone
params.Psk = 0.6265;% skin
params.Pkd = 2.1725;% kidney
params.Psp = 1.395;% spleen
params.Pgu = 1.0781;% gut
params.Pli = 1.9646;% liver
params.Pln = 1.2081;% lymph node
params.Poth = 1.0070;% others (median of all other values)


% Clearance and other constants
params.kF = 0.252; % kF = gut lumen transit rate (1/hr)
params.fR = 0.07; % for RIF
params.CL = 7.79; % systemic clearance rate RIF (L/h)  
params.ka = 1.07; %absorption rate RIF (1/hr)
params.kr = 0.17; %gut reabsorption rate for RIF (1/hr)

end

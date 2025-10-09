function all_params = loadPhysParams(drug, bw_PD, vol_PDs, vol_frac_PDs, qc_PD, flow_PDs, flow_frac_PDs)
% LOADPHYSPARAMS: Calculates random physiological parameter values based on
% given probability distributions
%
% DESCRIPTION:
%   This function repeatedly samples from a given array of probability
%   distributions in order to get a set of random parameters that describes 
%   a certain patient. It handles the distributions for raw values and
%   distributions that are for fractions of body weight or cardiac output
%   separately, making sure to renormalize the sampled fractions so their 
%   sum is equal to 1.
%
% INPUTS:
% - vol_PDs (cell array): Contains all probability distributions for raw
%   volume parameters not based on body weight
% - vol_frac_PDs (cell array): Contains all probability distributions for
%   volume parameters that are calculated as fractions of body weight 
% - flow_PDs (cell array): Contains all probability distributions for raw
%   flow parameters not based on cardiac output
% - flow_frac_PDs (cell array): Contains all probability distributions for
%   flow parameters that are calculated as fractions of cardiac output
%
% OUTPUTS:
% - all_params (cell array): Contains tables of compartment volume, blood
%   flow, lymph flow, and partition coefficient parameters (in that order)


%% Volumes

BW = random(bw_PD);

% sample all params from their probability dists
rand_params = zeros(1, length(vol_PDs));
for param_idx = 1:length(rand_params)
    rand_params(param_idx) = random(vol_PDs{param_idx});
end

rand_params_fracs = zeros(1, length(vol_frac_PDs));
for param_idx = 1:length(rand_params_fracs)
    rand_params_fracs(param_idx) = random(vol_frac_PDs{param_idx});
end

% renormalize fractions
vfracs_sum = sum(rand_params_fracs);
rand_params_fracs = (rand_params_fracs / vfracs_sum) * 0.870409878; % correction for handling some values raw
rand_params = [rand_params, rand_params_fracs .* BW];

% unpack and place into table
BW          = BW;
venblood    = rand_params(1);
artblood    = rand_params(2);
lymph       = rand_params(3);
lung        = rand_params(4);
brain       = rand_params(5);
heart       = rand_params(6);
adipose     = rand_params(7);
muscle      = rand_params(8);
skin        = rand_params(9);
kidney      = rand_params(10);
bone        = rand_params(11);
spleen      = rand_params(12);
gut         = rand_params(13);
liver       = rand_params(14);
other       = rand_params(15);
pleura      = rand_params(16);
fracsum     = vfracs_sum; % record sum of fractions for downstream analysis

vol_table = table(BW, venblood, artblood, lymph, ...
                    lung, brain, heart, adipose, muscle, ...
                    skin, kidney, bone, spleen, gut, ...
                    liver, other, pleura, fracsum);

%% Blood flows

QC = random(qc_PD);

% sample all params from their probability dists
rand_params = zeros(1, length(flow_PDs));
for param_idx = 1:length(rand_params)
    rand_params(param_idx) = random(flow_PDs{param_idx});
end

rand_params_fracs = zeros(1, length(flow_frac_PDs));
for param_idx = 1:length(rand_params_fracs)
    rand_params_fracs(param_idx) = random(flow_frac_PDs{param_idx});
end

% renormalize fractions
qfracs_sum = sum(rand_params_fracs);
rand_params_fracs = rand_params_fracs / qfracs_sum;
rand_params = [rand_params, rand_params_fracs .* QC];

% unpack and place into table
QC         = QC;
ka         = rand_params(1);
kdiss      = rand_params(2);
kF         = rand_params(3);
fR         = rand_params(4);
kr         = rand_params(5);
CL         = rand_params(6);
pleura     = rand_params(7);
bELF       = rand_params(8);
aELF       = rand_params(9);
artblood   = rand_params(10);
spleen     = rand_params(11);
gut        = rand_params(12);
brain      = rand_params(13);
heart      = rand_params(14);
adipose    = rand_params(15);
muscle     = rand_params(16);
bone       = rand_params(17);
skin       = rand_params(18);
other      = rand_params(19);
kidney     = rand_params(20);
liver      = spleen + gut + artblood;
fracsum    = qfracs_sum; % record sum of fractions for downstream analysis

flow_table = table(QC, ka, kdiss, fR, kF, kr, ...
                    CL, pleura, bELF, aELF, ...
                    artblood, spleen, gut, brain, ...
                    heart, adipose, muscle, bone, skin, ...
                    other, kidney, liver, fracsum);


%% Lymph flow in L/h
lymph     = 8/24; % L/day to L/h
lung      = lymph * 0.03;
brain     = lymph * 0.0105;
heart     = lymph * 0.01;
adipose   = lymph * 0.128;
muscle    = lymph * 0.16;
bone      = 0;
skin      = lymph * 0.0703;
other     = lymph * 0.0562;
kidney    = lymph * 0.085;
liver     = lymph * 0.33;
spleen    = 0;
gut       = lymph * 0.12;

lflow_table = table(lymph, lung, brain, heart, ...
                    adipose, muscle, bone, skin, ...
                    other, kidney, liver, spleen, ...
                    gut);

%% Partition coefficients (drug-specific)

if drug == "RIF"
    lung    = 1.7115;
    lymph   = 1.2081;
    brain   = 0.2285;
    adipose = 0.1885;
    heart   = 1.0158;
    muscle  = 0.6949;
    skin    = 0.6265;
    other   = 1.047;
    bone    = 0.3157;
    spleen  = 1.3950;
    kidney  = 2.1725;
    gut     = 1.0781;
    liver   = 1.9646;

elseif drug == "INHSA" || drug == "INHFA" % partition coeffs. are the same for either
    lung    = 0.7662;
    lymph   = 0.7556;
    brain   = 0.7537;
    adipose = 0.1543;
    heart   = 0.7551;
    muscle  = 0.7208;
    skin    = 0.6675;
    other   = 0.7435;
    bone    = 0.4330;
    spleen  = 0.7605;
    kidney  = 0.7441;
    gut     = 0.7429;
    liver   = 0.7212;

elseif drug == "PZA"
    lung    = 0.7381;
    lymph   = 0.7210;
    brain   = 0.7184;
    adipose = 0.1503;
    heart   = 0.7243;
    muscle  = 0.6868;
    skin    = 0.6496;
    other   = 0.7133;
    bone    = 0.4163;
    spleen  = 0.726;
    kidney  = 0.7127;
    gut     = 0.714;
    liver   = 0.6887;

elseif drug == "EMB"
    lung    = 4.3966;
    lymph   = 3.6743;
    brain   = 1.8054;
    adipose = 0.4625;
    heart   = 3.0623;
    muscle  = 2.7093;
    skin    = 1.9938;
    other   = 3.1337;
    bone    = 1.3765;
    spleen  = 4.0004;
    kidney  = 5.3375;
    gut     = 3.2051;
    liver   = 5.0703;

else
    error("Partition coeffs. not defined for the requested drug")

end

ptc_table = table(lung, lymph, brain, adipose, heart, ...
                    muscle, skin, other, bone, spleen, ...
                    kidney, gut, liver);

% package all param tables into cell array
all_params = {vol_table, flow_table, lflow_table, ptc_table};


end
function pt = loadPartitionCoefficients(drug)
% Load partition coefficients from Table S5
if strcmpi(drug, 'rifampicin')
    pt.Lu = 1.7115;
    pt.LN = 1.2081;
    pt.Tissue.Brain   = 0.2285;
    pt.Tissue.Adipose = 0.1885;
    pt.Tissue.Heart   = 1.0158;
    pt.Tissue.Muscle  = 0.6949;
    pt.Tissue.Skin    = 0.6265;
    pt.Tissue.Others  = 1.047;
    pt.Tissue.Bone    = 0.3157;
    pt.Tissue.Spleen  = 1.3950;
    pt.Tissue.Kidney  = 2.1725;
    pt.Tissue.Gut     = 1.0781;
    pt.Tissue.Liver   = 1.9646;
    
elseif strcmpi(drug, 'ethambutol')
    pt.Lu = 4.3966;
    pt.LN = 3.6743;
    pt.Tissue.Brain   = 1.8054;
    pt.Tissue.Adipose = 0.4625;
    pt.Tissue.Heart   = 3.0623;
    pt.Tissue.Muscle  = 2.7093;
    pt.Tissue.Skin    = 1.9938;
    pt.Tissue.Others  = 3.1337;
    pt.Tissue.Bone    = 1.3765;
    pt.Tissue.Spleen  = 4.0004;
    pt.Tissue.Kidney  = 5.3375;
    pt.Tissue.Gut     = 3.2051;
    pt.Tissue.Liver   = 5.0703;
    
 elseif strcmpi(drug, 'isoniazid')
    pt.Lu = 0.7662;
    pt.LN = 0.7556;
    pt.Tissue.Brain   = 0.7537;
    pt.Tissue.Adipose = 0.1543;
    pt.Tissue.Heart   = 0.7551;
    pt.Tissue.Muscle  = 0.7208;
    pt.Tissue.Skin    = 0.6675;
    pt.Tissue.Others  = 0.7435;
    pt.Tissue.Bone    = 0.4330;
    pt.Tissue.Spleen  = 0.7605;
    pt.Tissue.Kidney  = 0.7441;
    pt.Tissue.Gut     = 0.7429;
    pt.Tissue.Liver   = 0.7212;
    
elseif strcmpi(drug, 'pyrazinamide')
    pt.Lu = 0.7381;
    pt.LN = 0.7210;
    pt.Tissue.Brain   = 0.7184;
    pt.Tissue.Adipose = 0.1503;
    pt.Tissue.Heart   = 0.7243;
    pt.Tissue.Muscle  = 0.6868;
    pt.Tissue.Skin    = 0.6496;
    pt.Tissue.Others  = 0.7133;
    pt.Tissue.Bone    = 0.4163;
    pt.Tissue.Spleen  = 0.726;
    pt.Tissue.Kidney  = 0.7127;
    pt.Tissue.Gut     = 0.714;
    pt.Tissue.Liver   = 0.6887;
else
    error('Partition coefficients for %s not defined.', drug);
end
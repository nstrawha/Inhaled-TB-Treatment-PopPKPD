function KP = loadPartitionCoefficients_Lyons(drug)
% Partition coefficients from Lyons 2013 Table 2 
if strcmpi(drug, 'rifampin')
    KP.Lu = 0.49;
    KP.Ki  = 0.59;
    KP.Li  = 3.1;
    KP.Sp  = 0.25;
    KP.Gu  = 0.74;
    
    
    KP.Brain   = 0.11;
    KP.Adipose = 0.34;
    KP.Heart   = 0.74;
    KP.Muscle  = 0.76;
    KP.Skin    = 0.76;
    KP.Bone   = 0.22;
   
    KP.Oth = 0.59;
 
else
    error('Partition coefficients for %s not defined.', drug);
end
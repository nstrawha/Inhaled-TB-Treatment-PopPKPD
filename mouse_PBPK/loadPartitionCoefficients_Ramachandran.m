function KP = loadPartitionCoefficients_Ramachandran(drug)
% Load partition coefficients from Table S5
if strcmpi(drug, 'rifampicin')
%For RIF
KP.Lu = 1.7115; 
KP.Ki = 2.1725;
KP.Li = 1.9646;
KP.Sp = 1.3950;
KP.Gu = 1.0781;
KP.Oth = 1.047;

 
else
    error('Partition coefficients for %s not defined.', drug);
end
function KP = loadPartitionCoefficients(drug)
% These numbers are a mix of 
if strcmpi(drug, 'rifampicin')

%For RIF
KP.Lu = 0.441;   %Kplu = tissue to plasma partition coef for lung (Reali Supp Table 3) maybe we don't put this? It says the 'best fitted parameter'
KP.Ki = 2.1725; % Ramachandran Supp Table S5
KP.Li = 0.76; % where did this come from?
KP.Sp = 1.3950;% Ramachandran Supp Table S5
KP.Gu = 1.0781;% Ramachandran Supp Table S5

KP.Kpad  = 0.21; %adipose - where did this come from?
KP.Kpbr = 1; %brain - where did this come from?
KP.Kpbn = 0.3157; %bone  Ramachandran Supp Table S5
KP.Kpht = 1.0158; % heart Ramachandran Supp Table S5
KP.Kpm = 0.6949; % muscle Ramachandran Supp Table S5
KP.Kpsk = 0.6265;% skin Ramachandran Supp Table S5

KP.Oth = median(KP.Kpad + KP.Kpbr + KP.Kpbn + KP.Kpht + KP.Kpm  + KP.Kpsk); 
% AllKPs = [KP.Kpad; KP.Kpbr; KP.Kpbn; KP.Kpht; KP.Kpm; KP.Kpsk;KP.Lu;KP.Ki;KP.Li;KP.Sp;KP.Gu];
% KP.Oth = median(AllKPs);


else
    error('Partition coefficients for %s not defined.', drug);
end
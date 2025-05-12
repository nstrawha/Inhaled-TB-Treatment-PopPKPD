function phys = loadPhysiologyMouse(BW)

%% Blood flow rates
QC = 1.04; %(L/h) total cardiac output
% blood flow rate (fraction of total cardiac output) from Lee 2020 Supp
% Table 1
phys.Q.Sp = (0.011*QC); % spleen Lee 2020 - but should this include pancreas per Reali paper  
phys.Q.Ha = (0.02*QC); % hepatic artery Lee 2020
phys.Q.Li = (0.16*QC);% liver Lee 2020 matches Ruark 2014
phys.Q.Ki = (0.091*QC);%kidney Lee 2020
phys.Q.Gu = (0.11*QC); % gut but taking the intestines from  Lee 2020

%Other refers to 8 tissues: adipose, bone, brain, heart, muscle, gonads,
%and skin or is pancreas included here? Reali text is confusing for
%pancreas - didn't included gonads and pancreas for this though no numbers
% phys.QAd  = (0.070 *QC); %adipose  Lee 2020 matches fat Ruark 2014
% phys.QBr = (0.033*QC); %brain Lee 2020 matches Ruark 2014
% phys.QBn = 0.05*QC; %bone i think in Brown?
% phys.QHt = (0.066*QC); % heart Lee 2020
% phys.QMu = (0.16*QC); % muscle Lee 2020
% Qg = ; % gonads
% phys.QSk = 0.058;% skin
% phys.Qpa = (0.017*BW); %pancreas

% phys.Q.Oth = phys.QAd + phys.QBr + phys.QBn + phys.QHt + phys.QMu  + phys.QSk;
% OR can do the lungs minus what is known in this model: spleen, kidney,
% liver, gut 
phys.Q.Oth = 1- (0.011 + 0.16 + 0.091 + 0.11);

%% Tissue volumes - Lee 2020 mostly (Ruark and Brown for what's not in Lee 2020)
% BW is in kg
phys.V.A = 0.02;% arterial blood vol
phys.V.V = 0.045;% venous blood vol
phys.V.Sp = 0.0035*BW; %spleen 
phys.V.Gu = 0.042*BW;%gut
phys.V.Ki = 0.017*BW;%kidney
phys.V.Li = 0.055*BW;%liver
phys.V.Lu = 0.007*BW;%lung 

phys.Vad = 0.07*BW;
phys.Vbr = 0.017*BW;
phys.Vbn = 0.107*BW;
phys.Vht = 0.005*BW;
phys.Vmu = 0.384*BW;
phys.Vsk = 0.165*BW; 
% rest of body? 2.5 + 5.7 = 0.025 + 0.057

phys.V.Oth = phys.Vad + phys.Vbr + phys.Vbn + phys.Vht + phys.Vmu + phys.Vsk;

end


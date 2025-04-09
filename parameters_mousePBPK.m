function phys = loadPhysiologyMouse(BW)
%% RIF
phys.BP = 0.9;       %BP: blood to plasma ratio (Supp Table 1)

%% Blood flow rates
%Qlu,Qot,Qki,Qli;blood flow for: lung, lumped compart, kidneys and liver.
%Qsp, Qha,Qgu: blood flows for splenic, hepatic artery, and gut comp

QC = 1.04; %(L/h)
% blood flow rate (fraction of total cardiac output) from Lee 2020 Supp
% Table 1
phys.Q.Sp = (0.011*QC); % spleen - but should this include pancreas?  Lee 2020
phys.Q.Ha = (0.02*QC); % hepatic artery Lee 2020
phys.Q.Li = (0.16*QC);% liver Lee 2020 matches Ruark 2014
phys.Q.Ki = (0.091*QC);%kidney Lee 2020
phys.Q.Gu = (0.11*QC); % gut but taking the intestines?  Lee 2020

%Other refers to 8 tissues: adipose, bone, brain, heart, muscle, gonads,
%and skin or is pancreas included here? Reali text is confusing for
%pancreas - didn't included gonads and pancreasfor this though no numbers
phys.QAd  = (0.070 *QC); %adipose  Lee 2020 matches fat Ruark 2014
phys.QBr = (0.033*QC); %brain Lee 2020 matches Ruark 2014
phys.QBn = 0.05*QC; %bone i think in Brown?
phys.QHt = (0.066*QC); % heart Lee 2020
phys.QMu = (0.16*QC); % muscle Lee 2020
% Qg = ; % gonads
phys.QSk = 0.058;% skin
% phys.Qpa = (0.017*BW); %pancreas

phys.Q.Oth = phys.QAd + phys.QBr + phys.QBn + phys.QHt + phys.QMu  + phys.QSk;

%% Tissue volumes
% BW is in kg
phys.V.A = 0.02;% arterial blood vol
phys.V.V = .045;% venous blood vol
% total blood = .049
phys.V.Sp =0.0035*BW; %spleen 
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

phys.V.Oth = phys.Vad + phys.Vbr + phys.Vbn + phys.Vht + phys.Vmu + phys.Vsk ;


%% Partition coefficients - calcuated with Rowland Rodgers equations
% There are two papers with two equations.....Need
% ﻿the octanol:water partition coefficient (clogP), the vegetable oil:water partition coefficient (logD), the negative log10 of disassociation constants (pKa1, pKa2), affinity constant, the intracellular pH of red blood cells, the blood-to-plasma ratio (BP), and the free fraction in plasma (fup)—were extracted from the literature, internal sources, and DrugBank.

% clogP = 4.01; %Supp 1 this should be
% logD = 2.7; % Internet?
% pka1 = 1.7;%Supp 1
% pka2=7.9;%Supp 1
% affinity_constant = ;
% intrapH_RBCs = ;
% BP =BP;%Supp 1
% fup=fup; %Supp 1

% From Ramachandran 2023 who calculated using RR equations Table S5
%For RIF
phys.KP.Lu = 0.441;   %Kplu = tissue to plasma partition coef for lung (Supp Table 3) maybe we don't put this? It says the 'best fitted parameter'
% phys.KP.Lu = 1.7115;
phys.KP.Ki = 2.1725;
phys.KP.Li = 0.76;
phys.KP.Sp = 1.3950;
phys.KP.Gu = 1.0781;

phys.Kpad  = 0.21; %adipose
phys.Kpbr = 1; %brain
phys.Kpbn = 0.3157; %bone
phys.Kpht = 1.0158; % heart
phys.Kpm = 0.6949; % muscle
% Kpg = ; % gonads
phys.Kpsk = 0.6265;% skin
% Kppa = ; %pancreas

phys.KP.Oth = phys.Kpad + phys.Kpbr + phys.Kpbn + phys.Kpht + phys.Kpm  + phys.Kpsk;


%Kplu,Kpot,Kpki,Kpli;tissue to plasma partition coef for: lung, lumped compart, kidneys and liver
%Ksp,Kgu; partition coefficents for spleen and gut
end


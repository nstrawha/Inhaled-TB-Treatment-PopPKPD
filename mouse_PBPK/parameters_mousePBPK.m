function params = parameters_mousePBPK()
%% RIF
params.BP = 0.9;       %BP: blood to plasma ratio (Supp Table 1)
params.fup = 0.029;    %fup = fraction unbound in plasma (Supp Table 1)
params.CLr = 19;       %CLr = renal clearance fraction (Supp Table 1)
params.F = 1;          %F = drug bioavailability (Supp Table 1)
params.Ka = 0.3713;    %Ka = rate of absorption (Supp Table 3)It says the 'best fitted parameter'
params.CL = 0.037;     %Cl = total body clearance (Supp Table 3) It says the 'best fitted parameter'
params.Kplu = 0.441;   %Kplu = tissue to plasma partition coef for lung (Supp Table 3) maybe we don't put this? It says the 'best fitted parameter'

%% Blood flow rates
%Qlu,Qot,Qki,Qli;blood flow for: lung, lumped compart, kidneys and liver.
%Qsp, Qha,Qgu: blood flows for splenic, hepatic artery, and gut comp

% if body weight = 0.025 kg
% blood flow rate (fraction of total cardiac output) from Lee 2020 Supp
% Table 1
% where first term is the fraction of total cardiac output and 1.04 is cardiac output

params.Qsp = (0.011*1.04); % spleen - but should this include pancreas?  Lee 2020
params.Qha = (0.02*1.04); % hepatic artery Lee 2020
params.Qli = (0.16*1.04);% liver Lee 2020 matches Ruark 2014
params.Qki = (0.091*1.04);%kidney Lee 2020
params.Qgu = (0.11*1.04); % gut but taking the intestines?  Lee 2020

%Other refers to 8 tissues: adipose, bone, brain, heart, muscle, gonads,
%and skin or is pancreas included here? Reali text is confusing for
%pancreas
params.Qad  = (0.070 *1.04); %adipose  Lee 2020 matches fat Ruark 2014
params.Qbr = (0.033*1.04); %brain Lee 2020 matches Ruark 2014
params.Qbn = ; %bone i think in Brown?
params.Qht = (0.066*1.04); % heart Lee 2020
params.Qm = (0.16*1.04); % muscle Lee 2020
% Qg = ; % gonads
params.Qsk = ;% skin
params.Qpa = (0.017*1.04); %pancreas

params.Qot = params.Qad + params.Qbr + params.Qbn + params.Qht + params.Qm  + params.Qsk + params.Qpa;

%% Tissue volumes
params.Vab = ;%Vab = arterial blood vol
params.Vvb = ;%Vvb = venous blood vol
params.Vsp = %Vsp = spleen vol
params.Vgu = ;%Vgu = gut vol

params.Qot = params.Vad + params.Vbr + params.Vbn + params.Vht + params.Vm  + params.Vsk + params.Vpa;

params.Voth = 



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
params.Kplu = 1.7115;
params.Kpki = 2.1725;
params.Kpli = 0.76;
params.Ksp = 1.3950;
params.Kgu = 1.0781;

params.Kpad  = 0.21; %adipose
params.Kpbr = 1; %brain
params.Kpbn = 0.3157; %bone
params.Kpht = 1.0158; % heart
params.Kpm = 0.6949; % muscle
% Kpg = ; % gonads
params.Kpsk = 0.6265;% skin
% Kppa = ; %pancreas

params.Kpot = params.Kpad + params.Kpbr + params.Kpbn + params.Kpht + params.Kpm  + params.Kpsk;



%Kplu,Kpot,Kpki,Kpli;tissue to plasma partition coef for: lung, lumped compart, kidneys and liver
%Ksp,Kgu; partition coefficents for spleen and gut
end
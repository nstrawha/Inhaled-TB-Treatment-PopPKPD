
function Rif = RIF_Properties(B)
Rif.BP = 0.9; %Blood/plasma drug ratio
Rif.fu = 0.15; %Fraction unbound %%%%Have to find
Rif.KpuBC = round((B.H - 1 + Rif.BP)/(Rif.fu * B.H),2,"significant"); %Blood cell/plasma drug ratio
Rif.pKa1 = 1.7; %Acid dissociation constant 1
Rif.pKa2 = 7.9; %Acid dissociation constant 2
Rif.logPow = 2.7; %n-octanol:water partition coeff %%%%%%Humphries
Rif.logPvow = round(1.1115*Rif.logPow - 1.35,2,"significant"); %vegetable oil:water partition coeff
Rif.Pow = round(10^Rif.logPow,2,"significant"); %anti-log of logPow
Rif.Pvow = round(10^Rif.logPvow,2,"significant"); %anti-log of logPvow 

Rif.KaBC = round((Rif.KpuBC - (((1 + 10^(Rif.pKa2 - B.BCpH))/(1+10^(Rif.pKa2 - B.pHP)))*B.BCfracIW)...
- ((Rif.Pow*B.BCfracNL + (0.3*Rif.Pow + 0.7)*B.BCfracNP)/(1 + 10^(Rif.pKa2 - B.pHP))))...
*((1 + 10^(Rif.pKa2 - B.pHP))/(B.BCconcAP * 10^(Rif.pKa2 - B.BCpH))),2,"significant"); %Association constant of drug with acidic phospholipids in blood cells

%Pharmacokinetic parameters
Rif.CL = nan; %L/hr;
Rif.ka = nan; %1/hr %abs rate
Rif.fR = 0.183; %Fractional renal clearance
Rif.kF = 0.252; %1/hr %Gut lumen transit rate
Rif.kr = 0.17; %1/hr %Gut reabsorption rate

Rif.Name = "Rifampicin";
end
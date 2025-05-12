function Kp = calculateTissuePartition(Drug,TC,B)

% Tissue = {'Ad';'Bo';'Br';'Gu';'Hr';'Kd';'Li';'Lu';'Mu';'Sk';'Sp';'LN'};
% %original from Ramachandran model
Tissue = {'Ad';'Bo';'Br';'Gu';'Hr';'Ki';'Li';'Lu';'Mu';'Sk';'Sp';'LN'};

DrugP = [Drug.Pvow;Drug.Pow*ones(11,1)];

if Drug.Name=="Rifampicin"  %Zwitterion - 1
    X = 1 + 10^(Drug.pKa2 - TC.pHIW) + 10^(TC.pHIW - Drug.pKa1);
    Y = 1 + 10^(Drug.pKa2 - B.pHP) + 10^(B.pHP - Drug.pKa1);
    Kp_values = (TC.fracEW + (X.*TC.fracIW/Y) + (DrugP.*TC.fracNL + ((0.3*DrugP + 0.7).*TC.fracNP))/Y + ((Drug.KaBC*TC.concAP*(X-1))/Y))*Drug.fu; 

elseif Drug.Name=="Ethambutol"  %Moderate-to-strong base (Diprotic)
    X = 1 + 10^(Drug.pKa2 - TC.pHIW) + 10^(Drug.pKa1 + Drug.pKa2 - 2*TC.pHIW); 
    Y = 1 + 10^(Drug.pKa2 - B.pHP) + 10^(Drug.pKa1 + Drug.pKa2 - 2*B.pHP); 
    Kp_values = (TC.fracEW + (X.*TC.fracIW/Y) + (DrugP.*TC.fracNL + ((0.3*DrugP + 0.7).*TC.fracNP))/Y + ((Drug.KaBC*TC.concAP*(X-1))/Y))*Drug.fu; 

elseif Drug.Name=="Isoniazid"  %Weak base (Assumed to bind to albumin)
    X = 1 + 10^(Drug.pKa1 - TC.pHIW); 
    Y = 1 + 10^(Drug.pKa1 - B.pHP);
    Kp_values = (TC.fracEW + (X.*TC.fracIW/Y) + (DrugP.*TC.fracNL + ((0.3*DrugP + 0.7).*TC.fracNP))/Y... 
    + (((1/Drug.fu) - 1 - ((DrugP.*B.PlasmafracNL + (0.3*DrugP + 0.7).*B.PlasmafracNP)/Y)).*TC.RatioAlbPlasma))*Drug.fu;     
 
elseif Drug.Name=="Pyrazinamide"  %Neutral
    X = 1;
    Y = 1;
    Kp_values = (TC.fracEW + (X.*TC.fracIW/Y) + (DrugP.*TC.fracNL + ((0.3*DrugP + 0.7).*TC.fracNP))/Y... 
    + (((1/Drug.fu) - 1 - ((DrugP.*B.PlasmafracNL + (0.3*DrugP + 0.7).*B.PlasmafracNP)/Y)).*TC.RatioAlbPlasma))*Drug.fu; 

end

Kp = cell2struct(num2cell(Kp_values),Tissue,1);

Kp.Oth = median(Kp_values);
   
end

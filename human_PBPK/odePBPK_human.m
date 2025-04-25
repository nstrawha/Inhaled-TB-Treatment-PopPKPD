function dA = odePBPK(~, A, ka, kr,kF, CL, fR, phys, pt)
% Indexing convention for compartments:
% 1: Venous blood (CV), 2: Arterial blood (CA), 3: Lung (CLu)
% 4: Pleura (CPl), 5: Brain, 6: Adipose, 7: Heart, 8: Muscle
% 9: Skin, 10: Others, 11: Bone, 12: Spleen, 13: Kidney, 14: Gut
% 15: Liver, 16: Lymph Node, 17: Gut Lumen (AGL), 18: Absorbed drug (AD)

% Extract concentrations
CV = A(1); CA = A(2); CLu = A(3); CPl = A(4);
CT = A(5:15);  % Tissue compartments
CLN = A(16);   % Lymph node
AGL = A(17);   % Gut lumen
AD = A(18);    % Absorbed drug

% Parameters
Q = phys.Q;
L = phys.L;
PT_Tissue = struct2array(pt.Tissue);  % 11 tissue partition coefficients

% Calculate exiting concentrations from tissues
% CVT = CT ./ PT_Tissue(:);  % Ensure column vector
CVLu = CLu / pt.Lu;
CVLN = CLN / pt.LN;


dA = zeros(18,1);

%% Venous Blood
% Match CVT length with tissues in Q.Tissue and L.Tissue
% going to venous blood brain, adipose, heart, muscle, bone-0, skin, others,
% liver, kidney

CVT_V = [CT(1)/pt.Tissue.Brain,CT(2)/pt.Tissue.Adipose,CT(3)/pt.Tissue.Heart,CT(4)/pt.Tissue.Muscle,CT(5)/pt.Tissue.Skin,CT(6)/pt.Tissue.Others,CT(7)/pt.Tissue.Bone, CT(9)/pt.Tissue.Kidney,CT(11)/pt.Tissue.Liver];
        
dA(1) = sum((Q.TissuesV(1:9) - L.TissuesV(1:9)) .* CVT_V(1:9)) ...
    + (L.LN * CVLN) - (Q.total * CV);

%% Arterial Blood 
% arterial blood goes to brain, adipose, heart, muscle, bone, skin, others,
% gut,hepatic artery in place of liver, spleen, kidney 
% all 3 equations below give same results

%dA(2) = (Q.total-L.Lu) * CVLu - sum(Q.TissuesA) * CA; 
% dA(2) = (Q.total-L.Lu) * CVLu - Q.total * CA;
dA(2) = (Q.total-L.Lu) * CVLu - (Q.total-L.Lu) * CA;

%% Lung
dA(3) = Q.total * CV - (Q.total - L.Lu) * CVLu - (L.Lu - Q.Pl) * CVLu - Q.Pl * CVLu;

%% Pleura
dA(4) = Q.Pl * CVLu - Q.Pl * CPl;

%% Tissues with lymph flow (Brain, Adipose,Heart, Muscle, Skin, Others)
for i = 1:6
    idx = i + 4; % A(5) to A(10)
    Pi = PT_Tissue(i);
    CVTi = CT(i) / Pi;
    dA(idx) = Q.TissuesV(i) * CA - (Q.TissuesV(i) - L.TissuesV(i)) * CVTi - L.TissuesV(i) * CVTi;
end

%% Bone and Spleen (no lymph)

   
    CVTi_bo = CT(7) / pt.Tissue.Bone;
    CVTi_sp = CT(8) / pt.Tissue.Spleen;
    dA(11) = Q.Bo * CA - Q.Bo * CVTi_bo;
    dA(12) = Q.Sp * CA - Q.Sp * CVTi_sp;


%% Kidney (renal clearance)

i = 9; idx = 13;
Pi_kd = pt.Tissue.Kidney;
CVTi_kd = CT(i) / Pi_kd;
dA(idx) = Q.Kd * CA - (Q.Kd - L.Kd) * CVTi_kd - L.Kd * CVTi_kd - (fR * CL) * CA;

%% Gut (absorption + reabsorption)
i = 10; idx = 14;
Pi_gut = pt.Tissue.Gut;
CVT_Gu = CT(i) / Pi_gut;
dA(idx) = Q.Gu * CA - ((Q.Gu - L.Gu) * CVT_Gu) - L.Gu * CVT_Gu + ka * AD + kr * AGL;

%% Liver (hepatic clearance)
CVSp = CT(8) / pt.Tissue.Spleen;
CVGu = CT(10) / pt.Tissue.Gut;
CVLi = CT(11) / pt.Tissue.Liver;

QLi_in = (Q.LA * CA) + (Q.Sp * CVSp) + (Q.Gu-L.Gu) * CVGu;
QLi_out = (Q.Li - L.Li) * CVLi;

dA(15) = QLi_in - QLi_out - L.Li * CVLi - (1 - fR) * CL * (QLi_in / Q.Li);

%% Lymph Node
% going to lymph lung, brain, adipose, heart, muscle, skin, others,
% gut,liver, kidney
CVT_L=[CT(1)/pt.Tissue.Brain,CT(2)/pt.Tissue.Adipose, CT(3)/pt.Tissue.Heart,CT(4)/pt.Tissue.Muscle,CT(5)/pt.Tissue.Skin,CT(6)/pt.Tissue.Others,CT(9)/pt.Tissue.Kidney,CT(10)/pt.Tissue.Gut,CT(11)/pt.Tissue.Liver];
    
% CVT_L = CT([1,2,3,4,5,6,9,10,11]) ./ PT_Tissue([1,2,3,4,5,6,9,10,11]);


dA(16) = sum(L.TissuesL(1:9) .* CVT_L(1:9))+((L.Lu-Q.Pl)*CVLu)+(Q.Pl*CPl) - L.LN * CVLN;

%% Gut Lumen
dA(17) = (1 - fR) * CL * (QLi_in / Q.Li) - kr * AGL - kF * AGL;

%% Absorbed Drug
dA(18) = -ka * AD;
end

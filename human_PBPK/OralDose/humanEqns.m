function dA = humanEqns(~, A, ka, kr,kF, CL, fR, phys, pt)
% Indexing convention for compartments:
% 1: Venous blood, 2: Arterial blood, 3: Lung
% 4: Pleura, 5: Brain, 6: Adipose, 7: Heart, 8: Muscle
% 9: Skin, 10: Others, 11: Bone, 12: Spleen, 13: Kidney, 14: Gut
% 15: Liver, 16: Lymph Node, 17: Gut Lumen, 18: Absorbed drug

% Parameters
Q = phys.Q;
L = phys.L;
V = phys.V;

dA = zeros(18,1);

%% Venous Blood
dA(1) = (1/V.V)*((Q.Br - L.Br)*A(5)/pt.Tissue.Brain + (Q.Ad - L.Ad)*A(6)/pt.Tissue.Adipose + (Q.Hr - L.Hr) * A(7)/pt.Tissue.Heart + ...
    (Q.Mu-L.Mu) * A(8)/pt.Tissue.Muscle + (Q.Sk -L.Sk) * A(9)/pt.Tissue.Skin + (Q.Oth - L.Oth)*A(10)/pt.Tissue.Others + ...
    (Q.Bo)*A(11)/pt.Tissue.Bone + (Q.Kd-L.Kd)*A(13)/pt.Tissue.Kidney + (Q.Li-L.Li)*A(15)/pt.Tissue.Liver +...
    (L.LN * (A(16)/pt.LN)) - (Q.total * A(1)));

%% Arterial Blood
dA(2) = (1/V.A) * ((Q.total-L.Lu) * A(3)/pt.Lu - (Q.total-L.Lu) * A(2));

%% Lung
dA(3) = (1/V.Lu) * (Q.total * A(1) - (Q.total - L.Lu) * A(3)/pt.Lu - (L.Lu - Q.Pl) * A(3)/pt.Lu - Q.Pl * A(3)/pt.Lu);

%% Pleura
dA(4) = (1/V.Pl) * (Q.Pl * A(3)/pt.Lu - Q.Pl * A(4));

%% Tissues with lymph flow (Brain, Adipose,Heart, Muscle, Skin, Others)

dA(5) = (1/V.Brain) *(Q.Br * A(2) - (Q.Br - L.Br) * (A(5)/pt.Tissue.Brain) - L.Br * (A(5)/pt.Tissue.Brain));
dA(6) = (1/V.Adipose) *(Q.Ad * A(2) - (Q.Ad - L.Ad) * (A(6)/pt.Tissue.Adipose) - L.Ad * (A(6)/pt.Tissue.Adipose));
dA(7) = (1/V.Heart) *(Q.Hr * A(2) - (Q.Hr - L.Hr) * (A(7)/pt.Tissue.Heart) - L.Hr * (A(7)/pt.Tissue.Heart));
dA(8) = (1/V.Muscle) *(Q.Mu * A(2) - (Q.Mu - L.Mu) * (A(8)/pt.Tissue.Muscle) - L.Mu * (A(8)/pt.Tissue.Muscle));
dA(9) = (1/V.Skin) *(Q.Sk * A(2) - (Q.Sk - L.Sk) * (A(9)/pt.Tissue.Skin) - L.Sk * (A(9)/pt.Tissue.Skin));
dA(10) = (1/V.Others) *(Q.Oth * A(2) - (Q.Oth - L.Oth) * (A(10)/pt.Tissue.Others) - L.Oth * (A(10)/pt.Tissue.Others));

%% Bone and Spleen (no lymph)
dA(11) = (1/V.Bone) * (Q.Bo * A(2) - Q.Bo * (A(11)/pt.Tissue.Bone));
dA(12) = (1/V.Spleen) * (Q.Sp * A(2) - Q.Sp * (A(12)/pt.Tissue.Spleen));

%% Kidney (renal clearance)
dA(13) = (1/V.Kidney) * (Q.Kd * A(2) - (Q.Kd - L.Kd) *  A(13)/pt.Tissue.Kidney - L.Kd *  A(13)/pt.Tissue.Kidney - (fR * CL) * A(2));

%% Gut (absorption + reabsorption)
dA(14) = (1/V.Gut) * (Q.Gu * A(2) - (Q.Gu - L.Gu) * A(14)/pt.Tissue.Gut - L.Gu *  A(14)/pt.Tissue.Gut + ka * A(18) + kr * A(17));

%% Liver (hepatic clearance)

dA(15) = (1/V.Liver) * (Q.LA * A(2) + Q.Sp * A(12)/pt.Tissue.Spleen + (Q.Gu-L.Gu) * A(14)/pt.Tissue.Gut -  ...
    (Q.Li - L.Li) * A(15)/pt.Tissue.Liver - L.Li * A(15)/pt.Tissue.Liver - (1 - fR) * CL * ...
    (Q.LA * A(2) +Q.Sp * A(12)/pt.Tissue.Spleen + (Q.Gu-L.Gu) * A(14)/pt.Tissue.Gut) / Q.Li);

%% Lymph Node
dA(16) = (1/V.LN) * (L.Br * A(5)/pt.Tissue.Brain + L.Ad * A(6)/pt.Tissue.Adipose + L.Hr * A(7)/pt.Tissue.Heart +...
    L.Mu * A(8)/pt.Tissue.Muscle + L.Sk * A(9)/pt.Tissue.Skin + L.Oth * A(10)/pt.Tissue.Others+...
    L.Kd * A(13)/pt.Tissue.Kidney + L.Gu * A(14)/pt.Tissue.Gut + L.Li * A(15)/pt.Tissue.Liver + ...
    + (L.Lu-Q.Pl)*A(3)/pt.Lu +Q.Pl*A(4) - L.LN * A(16)/pt.LN);

%% Gut Lumen
dA(17) = (1 - fR) * CL * ((Q.LA * A(2) + Q.Sp * A(12)/pt.Tissue.Spleen + (Q.Gu-L.Gu) * A(14)/pt.Tissue.Gut) / Q.Li) - kr * A(17) - kF * A(17);

%% Absorbed Drug
dA(18) = -ka * A(18);
end

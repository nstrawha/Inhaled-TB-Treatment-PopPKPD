function dC = RIF_lung_ODEs(~, A, kdiss, kr, kF, effRB, effRA, br_frac, CL, fR, phys, pt)
% Script to set up human ODEs in all compartments

% Indexing convention for compartments:
% 1: Venous blood
% 2: Arterial blood
% 3: Lung
% 4: Pleura
% 5: Brain
% 6: Adipose
% 7: Heart
% 8: Muscle
% 9: Skin
% 10: Others
% 11: Bone
% 12: Spleen
% 13: Kidney
% 14: Gut
% 15: Liver
% 16: Lymph Node
% 17: Gut Lumen
% 18: Absorbed drug

%% Parameters

Q = phys.Q; % Tissue:blood flow rates
L = phys.L; % Tissue:lymph flow rates
V = phys.V; % Volumes of compartments

dC = zeros(20, 1);

%% Compartment  ----------------------------------------  Terms

% Venous Blood                                          
dC(1) = (1/V.V) * (...
        (Q.Br - L.Br) * A(5)/pt.Tissue.Brain + ...      % brain --> v. blood
        (Q.Ad - L.Ad) * A(6)/pt.Tissue.Adipose + ...    % adipose tissue --> v. blood
        (Q.Hr - L.Hr) * A(7)/pt.Tissue.Heart + ...      % heart --> v. blood
        (Q.Mu - L.Mu) * A(8)/pt.Tissue.Muscle + ...     % muscle --> v. blood
        (Q.Sk - L.Sk) * A(9)/pt.Tissue.Skin + ...       % skin --> v. blood
        (Q.Oth - L.Oth) * A(10)/pt.Tissue.Others + ...  % other tissue --> v. blood
        (Q.Bo) * A(11)/pt.Tissue.Bone + ...             % bone --> v. blood
        (Q.Kd - L.Kd) * A(13)/pt.Tissue.Kidney + ...    % kidney --> v. blood
        (Q.Li - L.Li) * A(15)/pt.Tissue.Liver + ...     % liver --> v. blood
        (L.LN * A(16)/pt.LN) - ...                      % lymph --> v. blood
        (Q.total * A(1)));                              % v. blood --> lungs

% Arterial Blood
dC(2) = (1/V.A) * (...
        (Q.total - L.Lu) * A(3)/pt.Lu - ... % lung --> a. blood
        (Q.total - L.Lu) * A(2));           % a. blood --> all other tissue

% Lung
dC(3) = (1/V.Lu) * (...
        Q.bELF * effRB * A(18) + ...        % bELF --> lung
        Q.aELF * effRA * A(19) - ...        % aELF --> lung
        Q.bELF * A(3) - ...                 % lung --> bELF
        Q.aELF * A(3) + ...                 % lung --> aELF
        Q.total * A(1) - ...                % v. blood --> lung
        (Q.total - L.Lu) * A(3)/pt.Lu - ... % lung --> a. blood
        (L.Lu - Q.Pl) * A(3)/pt.Lu - ...    % lung --> lymph
        Q.Pl * A(3)/pt.Lu);                 % lung --> pleura

% Pleura
dC(4) = (1/V.Pl) * (...
        Q.Pl * A(3)/pt.Lu - ... % lung --> pleura
        Q.Pl * A(4));           % pleura --> lung

% Tissues with lymph flow (Brain, Adipose, Heart, Muscle, Skin, Others)
% Brain
dC(5) = (1/V.Brain) * (...
        Q.Br * A(2) - ...                               % a. blood --> brain
        (Q.Br - L.Br) * A(5)/pt.Tissue.Brain - ...      % brain --> v. blood
        L.Br * A(5)/pt.Tissue.Brain);                   % brain --> lymph

% Adipose tissue
dC(6) = (1/V.Adipose) * (...
        Q.Ad * A(2) - ...                               % a. blood --> adipose
        (Q.Ad - L.Ad) * A(6)/pt.Tissue.Adipose - ...    % adipose --> v. blood
        L.Ad * A(6)/pt.Tissue.Adipose);                 % adipose --> lymph

% Heart
dC(7) = (1/V.Heart) * (...
        Q.Hr * A(2) - ...                               % a. blood --> heart
        (Q.Hr - L.Hr) * A(7)/pt.Tissue.Heart - ...      % heart --> v. blood
        L.Hr * A(7)/pt.Tissue.Heart);                   % heart --> lymph

% Muscle
dC(8) = (1/V.Muscle) * (...
        Q.Mu * A(2) - ...                               % a. blood --> muscle
        (Q.Mu - L.Mu) * A(8)/pt.Tissue.Muscle - ...     % muscle --> v. blood
        L.Mu * A(8)/pt.Tissue.Muscle);                  % muscle --> lymph

% Skin
dC(9) = (1/V.Skin) * (...
        Q.Sk * A(2) - ...                           % a. blood --> skin
        (Q.Sk - L.Sk) * A(9)/pt.Tissue.Skin - ...   % skin --> v. blood
        L.Sk * A(9)/pt.Tissue.Skin);                % skin --> lymph

% Other tissue
dC(10) = (1/V.Others) * (...
        Q.Oth * A(2) - ...                              % a. blood --> others
        (Q.Oth - L.Oth) * A(10)/pt.Tissue.Others - ...  % others --> v. blood
        L.Oth * A(10)/pt.Tissue.Others);                % others --> lymph

% Tissues with no lymph flow (Bone and Spleen)
% Bone 
dC(11) = (1/V.Bone) * (...
        Q.Bo * A(2) - ...               % a. blood --> bone
        Q.Bo * A(11)/pt.Tissue.Bone);   % bone --> v. blood

% Spleen
dC(12) = (1/V.Spleen) * (...
        Q.Sp * A(2) - ...                   % a. blood --> spleen
        Q.Sp * A(12)/pt.Tissue.Spleen);     % spleen --> liver

% Kidney (renal clearance)
dC(13) = (1/V.Kidney) * (...
        Q.Kd * A(2) - ...                               % a. blood --> kidney
        (Q.Kd - L.Kd) * A(13)/pt.Tissue.Kidney - ...    % kidney --> v. blood
        L.Kd *  A(13)/pt.Tissue.Kidney - ...            % kidney --> lymph
        (fR * CL) * A(2));                              % kidney clearance (urine)

% Gut (absorption + reabsorption)
dC(14) = (1/V.Gut) * (...
        Q.Gu * A(2) - ...                           % a. blood --> gut
        (Q.Gu - L.Gu) * A(14)/pt.Tissue.Gut - ...   % gut --> v. blood
        L.Gu *  A(14)/pt.Tissue.Gut + ...           % gut --> lymph
        kr * A(17));                                % gut lumen --> gut

% Liver (hepatic clearance)
dC(15) = (1/V.Liver) * (...
        Q.LA * A(2) + ...                               % a. blood --> liver
        Q.Sp * A(12)/pt.Tissue.Spleen + ...             % spleen --> liver
        (Q.Gu - L.Gu) * A(14)/pt.Tissue.Gut - ...       % gut --> liver
        (Q.Li - L.Li) * A(15)/pt.Tissue.Liver - ...     % liver --> v. blood
        L.Li * A(15)/pt.Tissue.Liver - ...              % liver --> lymph
        (1 - fR) * CL * (Q.LA * A(2) + ...              % liver --> gut lumen
        Q.Sp * A(12)/pt.Tissue.Spleen + ...             % liver --> gut lumen
        (Q.Gu - L.Gu) * A(14)/pt.Tissue.Gut) / Q.Li);   % liver --> gut lumen

% Lymph Node
dC(16) = (1/V.LN) * (...
        L.Br * A(5)/pt.Tissue.Brain + ...       % brain --> lymph
        L.Ad * A(6)/pt.Tissue.Adipose + ...     % adipose --> lymph
        L.Hr * A(7)/pt.Tissue.Heart + ...       % heart --> lymph
        L.Mu * A(8)/pt.Tissue.Muscle + ...      % muscle --> lymph
        L.Sk * A(9)/pt.Tissue.Skin + ...        % skin --> lymph
        L.Oth * A(10)/pt.Tissue.Others + ...    % others --> lymph
        L.Kd * A(13)/pt.Tissue.Kidney + ...     % kidney --> lymph
        L.Gu * A(14)/pt.Tissue.Gut + ...        % gut --> lymph
        L.Li * A(15)/pt.Tissue.Liver + ...      % liver --> lymph
        (L.Lu - Q.Pl) * A(3)/pt.Lu + ...        % lungs --> lymph
        Q.Pl * A(4) - ...                       % pleura --> lymph
        L.LN * A(16)/pt.LN);                    % lymph --> v. blood

% Gut Lumen
dC(17) = (1 - fR) * CL * (...
        (Q.LA * A(2) + ...                                  % a. blood --> gut lumen
        Q.Sp * A(12)/pt.Tissue.Spleen + ...                 % spleen --> gut lumen?
        (Q.Gu - L.Gu) * A(14)/pt.Tissue.Gut) / Q.Li) - ...  % gut lumen --> liver
        kr * A(17) - ...                                    % gut lumen --> gut (reabsorption)
        kF * A(17);                                         % gut lumen clearance (feces)

% Bronchi Epithelial Lining Fluid (bELF)
dC(18) = kdiss * A(20) * br_frac + ...  % dissolution --> bELF
        Q.bELF * A(3) - ...             % lung --> bELF
        Q.bELF * effRB * A(18);         % bELF --> lung

% Alveolar Epithelial Lining Fluid (aELF)
dC(19) = kdiss * A(20) * (1 - br_frac) + ...    % dissolution --> aELF
        Q.aELF * A(3) - ...                     % lung --> ELF
        Q.aELF * effRA * A(19);                 % aELF --> lung

% Drug Dissolution
dC(20) = -kdiss * A(20);   % dissolution --> ELF

end
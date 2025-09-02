function dC = oralODEs(~, C, params, effRB, effRA)
% ORALODES - Sets up the lung dose ODEs for rifampin.
%
% INPUTS:
% - C (double array): Contains the state vector of current
%   concentrations/amounts of drug in each compartment
% - params (cell array): A list of drug-specific parameters packaged in
%   tables
%
% OUTPUTS:
% dC (double array): Derivative of the concentration state vector

% Indexing for compartments:
% 1: Venous blood
% 2: Arterial blood
% 3: Lung
% 4: Pleura
% 5: brain
% 6: Adipose
% 7: Heart
% 8: Muscle
% 9: skin
% 10: other
% 11: Bone
% 12: Spleen
% 13: Kidney
% 14: Gut
% 15: Liver
% 16: Lymph Node
% 17: Gut Lumen
% 18: Bronchiolar Epithelial Lining Fluid (bELF)
% 19: Alveolar Epithelial Lining Fluid (aELF)
% 20: Drug Absorption/Dissolution


%% Unpack parameters

V   = params{1}; % volume parameter table
F   = params{2}; % blood flow parameter table
L   = params{3}; % lymph flor parameter table
K   = params{4}; % partition coefficient parameter table

dC = zeros(20, 1);


%% Compartment  ----------------------------------------  Terms

% Venous Blood                                          
dC(1) = (1/V.venblood) * (...
        (F.brain - L.brain) * C(5)/K.brain + ...            % brain --> v. blood
        (F.adipose - L.adipose) * C(6)/K.adipose + ...      % adipose tissue --> v. blood
        (F.heart - L.heart) * C(7)/K.heart + ...            % heart --> v. blood
        (F.muscle - L.muscle) * C(8)/K.muscle + ...         % muscle --> v. blood
        (F.skin - L.skin) * C(9)/K.skin + ...               % skin --> v. blood
        (F.other - L.other) * C(10)/K.other + ...           % other tissue --> v. blood
        (F.bone) * C(11)/K.bone + ...                       % bone --> v. blood
        (F.kidney - L.kidney) * C(13)/K.kidney + ...        % kidney --> v. blood
        (F.liver - L.liver) * C(15)/K.liver + ...           % liver --> v. blood
        (L.lymph * C(16)/K.lymph) - ...                     % lymph --> v. blood
        (F.QC * C(1)));                                     % v. blood --> lungs

% Arterial Blood
dC(2) = (1/V.artblood) * (...
        (F.QC - L.lung) * C(3)/K.lung - ... % lung --> a. blood
        (F.QC - L.lung) * C(2));            % a. blood --> all other tissue

% Lung
dC(3) = (1/V.lung) * (...
        F.QC * C(1) - ...                       % v. blood --> lung
        (F.QC - L.lung) * C(3)/K.lung - ...     % lung --> a. blood
        (L.lung - F.pleura) * C(3)/K.lung - ... % lung --> lymph
        F.pleura * C(3)/K.lung);                % lung --> pleura

%F.bELF * effRB * C(18) + ...            % bELF --> lung
%        F.aELF * effRA * C(19) - ...            % aELF --> lung
%        F.bELF * C(3) - ...                     % lung --> bELF
%        F.aELF * C(3) + ...                     % lung --> aELF

% Pleura
dC(4) = (1/V.pleura) * (...
        F.pleura * C(3)/K.lung - ... % lung --> pleura
        F.pleura * C(4));            % pleura --> lung

% Tissues with lymph flow (brain, adipose, heart, muscle, skin, other)
% Brain
dC(5) = (1/V.brain) * (...
        F.brain * C(2) - ...                       % a. blood --> brain
        (F.brain - L.brain) * C(5)/K.brain - ...   % brain --> v. blood
        L.brain * C(5)/K.brain);                   % brain --> lymph

% Adipose tissue
dC(6) = (1/V.adipose) * (...
        F.adipose * C(2) - ...                            % a. blood --> adipose
        (F.adipose - L.adipose) * C(6)/K.adipose - ...    % adipose --> v. blood
        L.adipose * C(6)/K.adipose);                      % adipose --> lymph

% Heart
dC(7) = (1/V.heart) * (...
        F.heart * C(2) - ...                       % a. blood --> heart
        (F.heart - L.heart) * C(7)/K.heart - ...   % heart --> v. blood
        L.heart * C(7)/K.heart);                   % heart --> lymph

% Muscle
dC(8) = (1/V.muscle) * (...
        F.muscle * C(2) - ...                       % a. blood --> muscle
        (F.muscle - L.muscle) * C(8)/K.muscle - ... % muscle --> v. blood
        L.muscle * C(8)/K.muscle);                  % muscle --> lymph

% Skin
dC(9) = (1/V.skin) * (...
        F.skin * C(2) - ...                   % a. blood --> skin
        (F.skin - L.skin) * C(9)/K.skin - ... % skin --> v. blood
        L.skin * C(9)/K.skin);                % skin --> lymph

% Other tissue
dC(10) = (1/V.other) * (...
        F.other * C(2) - ...                       % a. blood --> other
        (F.other - L.other) * C(10)/K.other - ...  % other --> v. blood
        L.other * C(10)/K.other);                  % other --> lymph

% Tissues with no lymph flow (Bone and Spleen)
% Bone 
dC(11) = (1/V.bone) * (...
        F.bone * C(2) - ...       % a. blood --> bone
        F.bone * C(11)/K.bone);   % bone --> v. blood

% Spleen
dC(12) = (1/V.spleen) * (...
        F.spleen * C(2) - ...           % a. blood --> spleen
        F.spleen * C(12)/K.spleen);     % spleen --> liver

% Kidney (renal clearance)
dC(13) = (1/V.kidney) * (...
        F.kidney * C(2) - ...                        % a. blood --> kidney
        (F.kidney - L.kidney) * C(13)/K.kidney - ... % kidney --> v. blood
        L.kidney *  C(13)/K.kidney - ...             % kidney --> lymph
        (F.fR * F.CL) * C(2));                       % kidney clearance (urine)

% Gut (absorption + reabsorption)
dC(14) = (1/V.gut) * (...
        F.gut * C(2) - ...                   % a. blood --> gut
        (F.gut - L.gut) * C(14)/K.gut - ...  % gut --> v. blood
        L.gut *  C(14)/K.gut + ...           % gut --> lymph
        F.ka * C(20) + ...                   % drug absorption
        F.kr * C(17));                       % gut lumen --> gut

% Liver (hepatic clearance)
dC(15) = (1/V.liver) * (...
        F.artblood * C(2) + ...                         % a. blood --> liver
        F.spleen * C(12)/K.spleen + ...                 % spleen --> liver
        (F.gut - L.gut) * C(14)/K.gut - ...             % gut --> liver
        (F.liver - L.liver) * C(15)/K.liver - ...       % liver --> v. blood
        L.liver * C(15)/K.liver - ...                   % liver --> lymph
        (1 - F.fR) * F.CL * (F.artblood * C(2) + ...    % liver --> gut lumen
        F.spleen * C(12)/K.spleen + ...                 % liver --> gut lumen
        (F.gut - L.gut) * C(14)/K.gut) / F.liver);      % liver --> gut lumen

% Lymph Node
dC(16) = (1/V.lymph) * (...
        L.brain * C(5)/K.brain + ...            % brain --> lymph
        L.adipose * C(6)/K.adipose + ...        % adipose --> lymph
        L.heart * C(7)/K.heart + ...            % heart --> lymph
        L.muscle * C(8)/K.muscle + ...          % muscle --> lymph
        L.skin * C(9)/K.skin + ...              % skin --> lymph
        L.other * C(10)/K.other + ...           % other --> lymph
        L.kidney * C(13)/K.kidney + ...         % kidney --> lymph
        L.gut * C(14)/K.gut + ...               % gut --> lymph
        L.liver * C(15)/K.liver + ...           % liver --> lymph
        (L.lung - F.pleura) * C(3)/K.lung + ... % lungs --> lymph
        F.pleura * C(4) - ...                   % pleura --> lymph
        L.lymph * C(16)/K.lymph);               % lymph --> v. blood

% Gut Lumen
dC(17) = (1 - F.fR) * F.CL * (...
        (F.artblood * C(2) + ...                         % a. blood --> gut lumen
        F.spleen * C(12)/K.spleen + ...                  % spleen --> gut lumen?
        (F.gut - L.gut) * C(14)/K.gut) / F.liver) - ...  % gut lumen --> liver
        F.kr * C(17) - ...                               % gut lumen --> gut (reabsorption)
        F.kF * C(17);                                    % gut lumen clearance (feces)

% Bronchi Epithelial Lining Fluid (bELF)
dC(18) = F.bELF * C(3) - ...               % lung --> bELF
         F.bELF * effRB * C(18);           % bELF --> lung

% Alveolar Epithelial Lining Fluid (aELF)
dC(19) = F.aELF * C(3) - ...                       % lung --> ELF
         F.aELF * effRA * C(19);                   % aELF --> lung

% Drug Absorption
dC(20) = -F.ka * C(20);   % absorption --> gut


end
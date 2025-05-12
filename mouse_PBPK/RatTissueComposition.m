function RTC = RatTissueComposition()
                  %[Adipose, Bone,   Brain,   Gut,    Heart,  Kidney, Liver,  Lung,   Muscle, Skin,   Spleen;  Lymph Node (mean of spleen and thymus)]
RTC.fracNL =         [0.853;  0.0174; 0.0391; 0.0375; 0.0135; 0.0121; 0.0135; 0.0215; 0.0100; 0.0603; 0.0071;  (0.0071+0.0168)/2]; %neutral lipid

RTC.fracNP =         [0.0016; 0.0016; 0.0015; 0.0124; 0.0106; 0.0240; 0.0238; 0.0123; 0.0072; 0.0044; 0.0107;  (0.0107+0.0092)/2];%neutral phospholipid

RTC.fracEW =         [0.135;  0.100;  0.162;  0.282;  0.320;  0.273;  0.161;  0.336;  0.118;  0.382;  0.207;   (0.207+0.150)/2];

RTC.fracIW =         [0.017;  0.346;  0.620;  0.475;  0.456;  0.483;  0.573;  0.446;  0.630;  0.291;  0.579;   (0.579+0.626)/2]; %Rodgers, 2005
 
RTC.concAP =         [0.40;   0.67;   0.40;   2.41;   2.25;   5.03;   4.56;   3.91;   1.53;   1.32;   3.18;    (3.18+2.30)/2];

RTC.RatioAlbPlasma = [0.049;  0.100;  0.048;  0.158;  0.157;  0.130;  0.086;  0.212;  0.064;  0.277;  0.097;   (0.097+0.075)/2];

RTC.pHIW = 7; %pH of intracellular tissue water
end



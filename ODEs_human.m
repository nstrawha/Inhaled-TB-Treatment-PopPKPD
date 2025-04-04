
function dydt = ODEs(t, y, params)

% % Unpack initial conditions from y
Cv = y(1);
Ca= y(2);
Clu= y(3);
Cpl= y(4);
Cbr= y(5);
Chr= y(6);
Cad= y(7);
Cmu= y(8);
Csk= y(9);
Coth= y(10);
Cbo= y(11);
Csp= y(12);
Ckd= y(13);
Cgu= y(14);
Cli= y(15);
Agl= y(16);
Cln= y(17);
Ad= y(18);

% Compute drug concentrations exiting each tissue
Cvln = Cln/params.Pln;
Cvlu = Clu/params.Plu;
Cvbr = Cbr/params.Pbr;
Cvhr = Chr/params.Phr;
Cvad = Cad/params.Pad;
Cvmu = Cmu/params.Pmu;
Cvsk = Csk/params.Psk;
Cvoth = Coth/params.Poth;
Cvbo = Cbo/params.Pbo;
Cvsp = Csp/params.Psp;
Cvkd = Ckd/params.Pkd;
Cvgu = Cgu/params.Pgu;
Cvli = Cli/params.Pli;

%% Equations
% params.Qt = flow rate to tissue organ 'T'
% params.Lt = lymph flow rate params.fRom tissue/organ 'T'
% Ca = drug conc in arterial blood
% Qpl = flow rate of the pleura
% params.CL = total systemic params.CLearance
% Ft = params.fRaction of total params.CLearnace apportioned to T (if any)
% Cvt = drug conc exiting T with Cvt = Ct/Pt where
% Pt = tissue:blood partition coefficient for T
% Summation of blood flow rates is for all tissues except lungs
% Amount of drug in tissue T is At = Ct*Vt wehre Vt = volume of T
% Ad = amount of drug input to the gut
% params.ka = oral absorption rate
% params.kr = rifampicin gut reabsorption rate during enteohepatic circulation


Cvt = Cvbr + Cvad + Cvhr + Cvmu + Cvbo + Cvsk + Cvkd + Cvsp + Cvgu + Cvli + Cvln + Cvoth;


% Eq 1 - Dynamics of drug concentration in venous blood
temp1 = (params.Qt-params.Lt)*Cvt;
Cv = (1/params.Vv) * (temp1 + (params.Lln *Cvln)-params.Qc*Cv);

% Eq 2 - Dynamics of drug concentration in arterial blood
temp2  = params.Qt*Ca;
Ca = (1/params.Va) * ((params.Qc-params.Llu) * Cvlu - temp2);

% Eq 3 - Drug concentration in lung compartment
Clu = (1/params.Vlu) * (params.Qc*Cv - (params.Qc-params.Llu)*Cvlu - (params.Llu-params.Qpl) *Cvlu - params.Qpl * Cvlu);

% Eq 4 - Drug concentration in pleura compartment
Cpl =  (1/params.Vpl) * (params.Qpl * Cvlu - params.Qpl*Cpl);

% Eq 5 -10 Drug concentration in non-eliminating tissues/organs with afferent
% lymph (brain, heart, adipose, musparams.CLe, skin, and others)
Cbr = (1/params.Vbr)* (params.Qbr * Ca - (params.Qbr-params.Lbr)*Cvbr-params.Lbr *Cvbr); %brain Eq 5

Chr =(1/params.Vhr)* (params.Qhr * Ca - (params.Qhr-params.Lhr)*Cvhr-params.Lhr *Cvhr); %heart Eq 6

Cad =(1/params.Vad)* (params.Qad * Ca - (params.Qad-params.Lad)*Cvad-params.Lad *Cvad); %adipose Eq 7

Cmu =(1/params.Vmu)* (params.Qmu * Ca - (params.Qmu-params.Lmu)*Cvmu-params.Lmu *Cvmu); %musparams.CLe Eq 8

Csk =(1/params.Vsk)* (params.Qsk * Ca - (params.Qsk-params.Lsk)*Cvsk-params.Lsk *Cvsk); %skin Eq 9

Coth =(1/params.Voth)* (params.Qoth * Ca - (params.Qoth-params.Loth)*Cvoth-params.Loth *Cvoth); %others Eq 10

% Eq 11 - 12 Drug concentration in non-eliminating tissues/organs without
% afferent lymph (bone, spleen)
Cbo = (1/params.Vbo) * (params.Qbo*Ca - params.Qbo*Cvbo); % bone Eq 12

Csp = (1/params.Vsp) * (params.Qsp*Ca - params.Qsp*Cvsp); % spleen Eq 12

% Eq 13 - Drug concentration in kidney compartment
Ckd = (1/params.Vkd) * (params.Qkd*Ca - (params.Qkd - params.Lkd)*Cvkd - params.Lkd *Cvkd-params.fR*params.CL*Ca);

% Eq 14 - Drug concentration in gut compartment
Cgu = (1/params.Vgu) * (params.Qgu*Ca - (params.Qgu-params.Lgu)*Cvgu - params.Lgu*Cvgu + params.ka *Ad + params.kr*Agl);

% Eq 15 - Drug concentration in liver compartment
temp3 = params.Qla*Ca + params.Qsp*Csp + params.Qgu*Cgu;
Cli = (1/params.Vli) * params.Qla*Ca + params.Qsp*Cvsp + (params.Qgu-params.Lgu)*Cvgu - (params.Qli-params.Lli)*Cvli - params.Lli * Cvli- (1-params.fR)*params.CL*(temp3/params.Qli);

% Eq 16 - Drug concentration in gut lumen compartment
Agl = (1-params.fR) * params.CL*(temp3/params.Qli) - params.kr*Agl - params.kF*Agl;

% Eq 17 - Drug concentration in lymph node compartment
temp4 = params.Lt*Cvt;
Cln = (1/params.Vln) *temp4 - params.Lln*Cvln;

% Eq 18 - Drug absorption
Ad = -params.ka *Ad;

dydt=[...
    Cv;...
    Ca;...
    Clu;...
    Cpl;...
    Cbr;...
    Chr;...
    Cad;...
    Cmu;...
    Csk;...
    Coth;...
    Cbo;...
    Csp;...
    Ckd;...
    Cgu;...
    Cli;...
    Agl;...
    Cln;...
    Ad;...
    ];


end


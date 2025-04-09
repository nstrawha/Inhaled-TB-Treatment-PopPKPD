
% Reali model

function dydt = ODEs(t, y, params)
%% Equations

% % Unpack initial conditions from y
Cab = y(1);
Cvb= y(2);
Clu= y(3);
Cot= y(4);
Csp= y(5);
Cgu= y(6);
Cor= y(7);
Cli= y(8);
Cki= y(9);

% Parameters
Q = phys.Q;
L = phys.L;
PT_Tissue = struct2array(pt.Tissue);  % 11 tissue partition coefficients

% Calculate exiting concentrations from tissues
CVT = CT ./ PT_Tissue(:);  % Ensure column vector
CVLu = CLu / pt.Lu;
CVLN = CLN / pt.LN;

dA = zeros(9,1);

% Eq 1 - Dynamics of drug concentration in arterial circulatory system
temp1 = Clu/(params.Kplu);
Qlu = params.Qsp + params.Qha + params.Qgu + params.Qki + params.Qot; %blood flow for: lung, lumped compart, kidneys and liver
dA(1) = (1/params.Vab)*(params.Qlu*temp1*params.BP-params.Qlu*Cab);

% Eq 2 - Dynamics of drug concentration in venous circulatory system
temp2 = params.Qot*(Cot/params.Kpot)*params.BP;
temp3 = params.Qki*(Cki/params.Kpki)*params.BP;
temp4 = params.Qli*(Cli/params.Kpli)*params.BP;
dA(2)= (1/params.Vvb)*(temp2+temp3+temp4-params.Qlu*Cvb);

% Eq 3 - Drug concentration in lung compartment
temp5 = params.Qlu*(Clu/params.Kplu)*params.BP;
dA(3) = (1/params.Vlu)*(params.Qlu*Cvb-temp5);

% Eq 4 - Drug concentration in 'other' compartment
temp6 = (Cot/params.Kpot)*params.BP;
dA(4)= (1/params.Vot)*(params.Qot*(Cab-temp6));

% Eq 5 - Drug concentration in spleen compartment
temp7 = (Csp/params.Kpsp)*params.BP;
dA(5)=(1/params.Vsp)*(params.Qsp*(Cab-temp7));

% Eq 6 - Drug concentration in gut compartment
temp8 = params.Ka*Cor*params.F;
temp9 = (Cgu/params.Kpgu)*params.BP;
dA(6) = (1/params.Vgu)*(temp8 + params.Qgu * (Cab-temp9));

% Eq 7 - Drug concentration in oral dose compartment
dA(7) = -(params.Ka*Cor*params.F);

% Eq 8 - Drug concentration in liver compartment
temp10=params.Qha*Cab;
temp11=params.Qgu*(Cgu/params.Kpgu)*params.BP;
temp12=params.Qsp*(Csp/params.Kpsp)*params.BP;
temp13=params.Qli*(Cli/params.Kpli)*params.BP;
dA(8)=(1/params.Vli)*(temp10 + temp11 + temp12 - temp13 - params.CL * (1-params.CLr) * Cli * params.fup);

% Eq 9 - Drug concentration in kidney compartment
temp14=(Cki/params.Kpki)*params.BP;
dA(9) = (1/params.Vki) *(params.Qki * (Cab- temp14) - params.CL * params.Clr * Cki * params.fup);


end
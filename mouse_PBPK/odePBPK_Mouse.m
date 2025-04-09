
% Reali model

function dA = odePBPK_Mouse(~, y,fup,CLr,F,ka,CL, phys)
%% Equations

% % Unpack initial conditions from y
Cab = y(1); % arterial blood compartment
Cvb= y(2); % venous blood compartment
Clu= y(3); % lung compartment
Cot= y(4); % other compartment
Csp= y(5); % spleen compartment
Cgu= y(6); % gut compartment
Cor= y(7); % oral dose compartment
Cli= y(8); % liver compartment
Cki= y(9); % kidney compartment

% Parameters
BP = phys.BP;
Q = phys.Q;
V = phys.V;
KP = phys.KP;

dA = zeros(9,1);

% Eq 1 - Dynamics of drug concentration in arterial circulatory system
temp1 = Clu/(KP.Lu);
Qlu = Q.Sp + Q.Ha + Q.Gu + Q.Ki + Q.Oth; %blood flow for: lung, lumped compart, kidneys and liver
dA(1) = (1/V.A)*(Qlu*temp1*BP-Qlu*Cab);

% Eq 2 - Dynamics of drug concentration in venous circulatory system
temp2 = Q.Oth*(Cot/KP.Oth)*BP;
temp3 = Q.Ki*(Cki/KP.Ki)*BP;
temp4 = Q.Li*(Cli/KP.Li)*BP;
dA(2)= (1/V.V)*(temp2+temp3+temp4-Qlu*Cvb);

% Eq 3 - Drug concentration in lung compartment
temp5 = Qlu*(Clu/KP.Lu)*BP;
dA(3) = (1/V.Lu)*(Qlu*Cvb-temp5);

% Eq 4 - Drug concentration in 'other' compartment
temp6 = (Cot/KP.Oth)*BP;
dA(4)= (1/V.Oth)*(Q.Oth*(Cab-temp6));

% Eq 5 - Drug concentration in spleen compartment
temp7 = (Csp/KP.Sp)*BP;
dA(5)=(1/V.Sp)*(Q.Sp*(Cab-temp7));

% Eq 6 - Drug concentration in gut compartment
temp8 = ka*Cor*F;
temp9 = (Cgu/KP.Gu)*BP;
dA(6) = (1/V.Gu)*(temp8 + Q.Gu * (Cab-temp9));

% Eq 7 - Drug concentration in oral dose compartment
dA(7) = -(ka*Cor*F);

% Eq 8 - Drug concentration in liver compartment
temp10=Q.Ha*Cab;
temp11=Q.Gu*(Cgu/KP.Gu)*BP;
temp12=Q.Sp*(Csp/KP.Sp)*BP;
temp13=Q.Li*(Cli/KP.Li)*BP;
dA(8)=(1/V.Li)*(temp10 + temp11 + temp12 - temp13 - CL * (1-CLr) * Cli * fup);

% Eq 9 - Drug concentration in kidney compartment
temp14=(Cki/KP.Ki)*BP;
dA(9) = (1/V.Ki) *(Q.Ki * (Cab- temp14) - CL * CLr * Cki * fup);


end
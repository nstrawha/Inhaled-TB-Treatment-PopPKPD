
% Reali model

function dA = mouseEqns(~, A,fup,CLr,F,ka,CL,BP, phys, KP)
%% Equations

% % Unpack initial conditions from y
% 1 arterial blood compartment
% 2 venous blood compartment
% 3 lung compartment
% 4 other compartment
% 5 spleen compartment
% 6 gut compartment
% 7 oral dose compartment
% 8 liver compartment
% 9 kidney compartment

% Parameters
Q = phys.Q;
V = phys.V;

dA = zeros(9,1);

% Eq 1 - Dynamics of drug concentration in arterial circulatory system
dA(1) = (1/V.A)*((Q.Sp + Q.Ha + Q.Gu + Q.Ki + Q.Oth) * (A(3)/KP.Lu) * BP - (Q.Sp + Q.Ha + Q.Gu + Q.Ki + Q.Oth)*A(1));

% Eq 2 - Dynamics of drug concentration in venous circulatory system
dA(2)= (1/V.V)*( Q.Oth*(A(4)/KP.Oth)*BP + Q.Ki * (A(9)/KP.Ki) * BP +  Q.Li*(A(8)/KP.Li)*BP - (Q.Sp + Q.Ha + Q.Gu + Q.Ki + Q.Oth) * A(2));

% Eq 3 - Drug concentration in lung compartment
dA(3) = (1/V.Lu)*((Q.Sp + Q.Ha + Q.Gu + Q.Ki + Q.Oth)*A(2)-(Q.Sp + Q.Ha + Q.Gu + Q.Ki + Q.Oth)*(A(3)/KP.Lu)*BP);

% Eq 4 - Drug concentration in 'other' compartment
dA(4)= (1/V.Oth)*(Q.Oth*(A(1)-(A(4)/KP.Oth)*BP));

% Eq 5 - Drug concentration in spleen compartment
dA(5)=(1/V.Sp)*(Q.Sp*(A(1)-(A(5)/KP.Sp)*BP));

% Eq 6 - Drug concentration in gut compartment
dA(6) = (1/V.Gu)*(ka*A(7)*F + Q.Gu * (A(1)-(A(6)/KP.Gu)*BP));

% Eq 7 - Drug concentration in oral dose compartment
dA(7) = -(ka*A(7)*F);

% Eq 8 - Drug concentration in liver compartment
dA(8)=(1/V.Li)*(Q.Ha*A(1) + Q.Gu*(A(6)/KP.Gu)*BP + Q.Sp*(A(5)/KP.Sp)*BP - Q.Li*(A(8)/KP.Li)*BP - CL * (1-CLr) * A(8) * fup);

% Eq 9 - Drug concentration in kidney compartment
dA(9) = (1/V.Ki) *(Q.Ki * (A(1)- (A(9)/KP.Ki)*BP) - CL * CLr * A(9) * fup);


end
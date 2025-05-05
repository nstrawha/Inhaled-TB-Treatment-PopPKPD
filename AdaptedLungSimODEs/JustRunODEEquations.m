
% Read in parameters from CSV file in this folder
ODEParamsRead = readtable('ODEParams.csv');

% Change the parameters into a structure version instead of table and get
% the first (minimum) number
struct_versionMin = struct();
for i = 1:height(ODEParamsRead)
    struct_versionMin.(ODEParamsRead.ParameterName{i}) = ODEParamsRead.Min(i);
end

% Change the parameters into a structure version instead of table and get
% the first (maximum) number
struct_versionMax = struct();
for i = 1:height(ODEParamsRead)
    struct_versionMax.(ODEParamsRead.ParameterName{i}) = ODEParamsRead.Max(i);
end

% implement lhs here


t=0;
tEnd = 200;
tStart=t;
opts = odeset('NonNegative',1);
% opts is likely for non-negative (said by Christian)

ODEParams = struct_versionMin;

% start with 1 infected mac, 1 intracellular mtb and 1 activated mac
% can start with 0 everything else though
y = zeros(1,16); y(2) = 1; y(3) =1; y(11)=1;

[tTemp,yTemp]=ode15s(@(t,y) ODEquations(t,y,ODEParams),[tStart tEnd],y,opts);

t=tTemp(end);
yTemp=real(yTemp);
y=yTemp(end,:);

names= {'MR';...
    'MI';...
    'MA';...
    'T0';...
    'T1';...
    'T2';...
    'T80';...
    'TC';...
    'T8';...
    'BE';...
    'BI';...
    'IG';...
    'I12';...
    'I10';...
    'I4';...
    'TNF'};
figure;
tiledlayout(4,4);
set(gcf,'Position',[00 00 1920 1080])

for i=1:16
    nexttile
    plot(tTemp, yTemp(:,i),'LineWidth',2);
    title(names{i});
    hold on;
end

% Fig 4 
CFU = yTemp(:,10) + yTemp(:,11);
Macs = yTemp(:,1) + yTemp(:,2)+ yTemp(:,3);
CD3Tcells = yTemp(:,4) + yTemp(:,5)+ yTemp(:,6) + yTemp(:,7) + yTemp(:,8)+ yTemp(:,9);

figure;
% A Total CFU
set(gcf,'Position',[00 00 1920 1080])
tiledlayout(1,3)
nexttile
plot(tTemp, CFU,'LineWidth',2); 
set(gca, 'YScale', 'log');
title('Total CFU');
hold on;

% B Total Macrophages
nexttile
plot(tTemp, Macs,'LineWidth',2);
set(gca, 'YScale', 'log');
title('Total Macrophages');
hold on;

% C Total CD3 T cells
nexttile
plot(tTemp, CD3Tcells,'LineWidth',2);
set(gca, 'YScale', 'log');
title('Total CD3 T cells');
hold on;




% if ismember(t,StateStoreTimeList) %put this loop in function?????
%     Func_StoreCurrentState
% end

%{
description:
this is main algorithm of the code--it runs the ODEs

used by:
RunLHS.m

uses:
Func_RungeKutta4thOrder.m
ODEs_DeniseOriginal.m
ODEquations.m
Func_AddToPlotList.m
Func_RunPassFailTests.m


NOTE: should do something to make sure RK4 interval hits target evaluation
times (stepping by dt, and it's possible end time is between dt's)


NOTE: should rename and perhaps combine Func_AddMR1stIter.m and
Func_ComputeNewRestingMacs.m
in those, PropCompartmentsOccupied needs to be defined in Inputs, and
probably needs to be adjusted, and, in fact, should just have VolConc found
rather than computing using GranSim definitions.
also, determine more certainly whether or not methods are equivalent, and
if they are then remove the iteration method.
%}


ODESolveStartTime=toc;

TimeIndex=1;
IterCount=0;

t=StartTime;
y=y0;

%FOLLOWING REMOVED ONCE ADDED Func_StoreCurrentState oct 22 2018
% PlotIter=PlotIter+1;
% tPlot(PlotIter)=t;
% yPlot(PlotIter,:)=y;


OutputUpdateCounter=-1;


ContinueTest=0;
while (t<StopTime) && (ContinueTest==0)
    
    if TimeIndex<length(EvaluationTimesList)
        TimeIndex=TimeIndex+1;
        tEnd=EvaluationTimesList(TimeIndex);
        if tEnd>StopTime
            tEnd=StopTime;
        end
    else
        tEnd=StopTime;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %RK4 method
    %----------------------------------------------------------------------
    if (strcmp(ODEMethod,'RK4thOrder')==1)
        
        if strcmp(MR_GranSizeMethodYesOrNo,'yes')
            Func_AddMR1stIter
        end
        
        
        while (tEnd-t)>dt/2 %since adding floating pts should do this to prevent over-shooting
            
            
            OutputUpdateCounter=OutputUpdateCounter+1;
            if OutputUpdateCounter>CommandOutputIterStride-1
                OutputUpdateCounter=0;
                fprintf('\n%f\n',t)
            end
            
            
            [tNew,yNew]=Func_RungeKutta4thOrder(t,y,ODEParams);
            %==============================================================
            %try to make non-negative
            %--------------------------------------------------------------
            %             if min(yNew)<0
            %                 tStartTemp=t;tEndTemp=tNew;
            %                 ODEParams.ODEMethod='Stiff';
            %                 if strcmp(ModelDescription,'Denise')==1
            %                     [tTemp,yTemp]=ode15s(@(t,y) ODEs_DeniseOriginal(t,y,ODEParams),[tStartTemp tEndTemp],y,opts);
            %                 elseif strcmp(ModelDescription,'New')==1
            %                     %[tTemp,yTemp]=ode15s(@(t,y) ODEquations(t,y,ODEParams),[tStartTemp tEndTemp],y,opts);
            %                     [tTemp,yTemp]=ode15s(@(t,y) ODEquations_UnpackEmbedded(t,y,ODEParams),[tStartTemp tEndTemp],y,opts);
            %                 end
            %                 ODEParams.ODEMethod='RK4thOrder';
            %                 tNew=tTemp(end);
            %                 yNew=yTemp(end,:);
            %             end
            yNew=max(real(yNew),0);
            
            if strcmp(MR_GranSizeMethodYesOrNo,'yes')
                Func_ComputeNewRestingMacs
            end
            
            %==============================================================
            %Func_AddToPlotList %REMOVED WHEN ADDED Func_StoreCurrentState
            t=tNew;
            y=yNew;
        end
        t=tEnd; %since adding floating pts should do this to ensure hitting the exact number for testing
        if ismember(t,StateStoreTimeList)
            Func_StoreCurrentState
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %stiff method
    %----------------------------------------------------------------------
    if (strcmp(ODEMethod,'Stiff')==1) || (strcmp(ODEMethod,'Stiff_Events')==1)
        tStart=t;
        if strcmp(ModelDescription,'Denise')==1
            [tTemp,yTemp]=ode15s(@(t,y) ODEs_DeniseOriginal(t,y,ODEParams),[tStart tEnd],y,opts);
        elseif strcmp(ModelDescription,'New')==1
            [tTemp,yTemp]=ode15s(@(t,y) ODEquations(t,y,ODEParams),[tStart tEnd],y,opts);
            %[tTemp,yTemp]=ode15s(@(t,y) ODEquations_UnpackEmbedded(t,y,ODEParams),[tStart tEnd],y,opts);
%             [tTemp,yTemp]=ode23t(@(t,y) ODEquations_UnpackEmbedded(t,y,ODEParams),[tStart tEnd],y,opts);
            %yTemp=max(yTemp,0);
        end
        t=tTemp(end);
        yTemp=real(yTemp);
        y=yTemp(end,:);
        %removed below lines when added Func_StoreCurrentState
%         tPlot=[tPlot;tTemp(2:end)]; %since variable step size can't preallocate, but i think it's worth inefficiency
%         yPlot=[yPlot;yTemp(2:end,:)];
        if ismember(t,StateStoreTimeList) %put this loop in function?????
            Func_StoreCurrentState
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %stiff method PROPOSED CHANGED VERSION
    %----------------------------------------------------------------------
    if (strcmp(ODEMethod,'Stiff_NEW')==1)
        while t<tEnd
            tStart=t;
            tEndTemp=min(tStart+dtStiffInterval,tEnd);
            StartODETimer=toc;
            if strcmp(ModelDescription,'Denise')==1
                [tTemp,yTemp]=ode15s(@(t,y) ODEs_DeniseOriginal(t,y,ODEParams),[tStart tEndTemp],y,opts);
            elseif strcmp(ModelDescription,'New')==1
                [tTemp,yTemp]=ode15s(@(t,y) ODEquations(t,y,ODEParams),[tStart tEnd],y,opts);
%                 [tTemp,yTemp]=ode15s(@(t,y) ODEquations_UnpackEmbedded(t,y,ODEParams),[tStart tEnd],y,opts);
                %yTemp=max(yTemp,0);
            end
            EndODETimer=toc;
            if (EndODETimer-StartODETimer)>ODERunTimeThreshold
                PassFailODESolverClockTest='Fail';
                break
            end
            t=tTemp(end);
            yTemp=real(yTemp);
            y=yTemp(end,:);
            y=real(y);
        %removed below lines when added Func_StoreCurrentState
%         tPlot=[tPlot;tTemp(2:end)]; %since variable step size can't preallocate, but i think it's worth inefficiency
%         yPlot=[yPlot;yTemp(2:end,:)];
        if ismember(t,StateStoreTimeList) %put this loop in function?????
            Func_StoreCurrentState
        end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Func_RunPassFailTests
    
    
end


%fprintf('\nN=%f, BI/MI=%f',N,y(end,11)/y(end,2))
ODESolveRuntime=toc-ODESolveStartTime;

%RunTimeSampleList(SampleNumber)=ODESolveRuntime;
fprintf('\nFinished sample number %i\n Runtime: %f',SampleNumber,ODESolveRuntime)







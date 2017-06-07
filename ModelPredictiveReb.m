function [rebalanceQueue]=ModelPredictiveReb(RoadNetwork,RebWeight,Passengers,Flags)

% Taken from AMoD-power:TIBalancedPowerFlow and edited as needed.

% Function to solve time-invariant AMoD routing and rebalancing with
% charging

% Syntax: [cplex_out]=SolveOscarrBattery(RoadGraph,RoadCap,TravelTimes,ChargeReqs,ChargerTimes,C,RebWeight,Sources,Sinks,FlowsIn,PowerGraph,PowerLineCap,PowerLineResistance,PowerGensList,PowerGensMax,PowerGensMin,PowerCosts,PowerExtLoads,RoadToPowerMap,PowerToRoadMap,milpflag,congrelaxflag,sourcerelaxflag,cachedAeqflag)

% Inputs:
% * Inside RoadNetwork
% - RoadGraph, a nx1 cell structure. RoadGraph{i} contains the neighbors of
%    i in the road graph
% - RoadCap, a nxn matrix. RoadCap(i,j) is the capacity of the i-j link
%    (in vehicles per unit time).
% - TravelTimes, a nxn matrix. TravelTimes(i,j) is the travel time along
%    the i-j link.
% - ChargeReqs, a nxn matrix. ChargeReqs(i,j) is the amount of (quantized)
%    units of charge required to travel from i to j.
% - ChargerTimes, a nx1 matrix. ChargerTimes(i) is the time required to
%    charge by one unit charge at node i.
% - C, a scalar. The overall number of charge levels.
% - ChargerCaps, nx1 vector. ChargerCaps(i) returns the capacity of the charger at node i
% * RebWeight, the relative importance of the rebalancing cost wrt the
%    customer cost
% * Inside Passengers
% - Sources, a S-x-1 vector. Sources(i) is the ith source node
% - Sinks, Sx1 cell structure. Sinks{i} contains the sinks of the classes that begin at Source(i)
% - FlowsIn, a m-by-1 cell array with flows of class i. FlowsIn{i}(j) is the amount of (vehicle) flow
% entering the source i in class j.
% * Inside Flags
% - milpflag (default: 1). If 1, the vehicle routing problem is solved as a
%    MILP. If 0, the problem is solved as a linear relaxation.
% - congrelaxflag (default: 0). If this flag is on, then vehicle flows are
%    allowed to violate the  congestion constraint for a cost.
%    A slack variable is defined for each edge. The cost is defined in the
%    main code.
% - sourcerelaxflag (default: 0). %If this flag is on, each vehicle flow is
%    allowed to reduce its sources (and sinks) for a cost. This is
%    especially useful when it is not possible to compute a satisfying
%    rebalancing flow because of timing constraints but, if the congestion
%    constraints are tight, it can also preserve feasibility if not all
%    flows can be realized. A slack variable is defined for each source and sink.

%% Unpack things

T=RoadNetwork.T;
RoadGraph=RoadNetwork.RoadGraph;
TravelTimes=RoadNetwork.TravelTimes;
Starters=RoadNetwork.Starters; % Starters(i) gives the number of available vehs at t=1

FlowsIn=Passengers.FlowsIn;
FlowsOut = Passengers.FlowsOut;

milpflag=Flags.milpflag;
congrelaxflag=Flags.congrelaxflag;
sourcerelaxflag=Flags.sourcerelaxflag;
cachedAeqflag=Flags.cachedAeqflag;

%% Utilities

debugflag=1;        %Makes output verbose
DIAGNOSTIC_FLAG=0;  % Diagnoses state allocation. Useful for initial debugging.

%CongestionCost=1e3; %The cost of violating the congestion constraint
SourceRelaxCost=1e6;     % The cost of dropping a source or sink altogether

%Clean up road graph.
for i=1:length(RoadGraph)
    RoadGraph{i}=sort(unique(RoadGraph{i}));
end

%Nodes in ReverseRoadGraph{i} are such that RoadGraph{ReverseRoadGraph{i}} contains
%i
ReverseRoadGraph=cell(size(RoadGraph));
for i=1:length(RoadGraph)
    for j=RoadGraph{i}
        ReverseRoadGraph{j}=[ReverseRoadGraph{j} i];
    end
end
for i=1:length(ReverseRoadGraph)
    ReverseRoadGraph{i}=sort(unique(ReverseRoadGraph{i}));
end
%Nodes in ReversePowerGraph{i} are such that PowerGraph{ReversePowerGraph{i}} contains
%i


N=length(RoadGraph);
M=1;
S = N;

E=0;
NumRoadEdges=zeros(N,1);
for i=1:N
    NumRoadEdges(i)=length(RoadGraph{i});
    E=E+length(RoadGraph{i});
end
cumRoadNeighbors=cumsum(NumRoadEdges);
cumRoadNeighbors=[0;cumRoadNeighbors(1:end-1)];

RoadNeighborCounter=sparse([],[],[],N,N,E);
TempNeighVec=zeros(N,1);
for i=1:N
    for j=RoadGraph{i}
        TempNeighVec(j)=1;
    end
    NeighCounterLine=cumsum(TempNeighVec);
    for j=RoadGraph{i}
        RoadNeighborCounter(i,j)=NeighCounterLine(j);
    end
    TempNeighVec=zeros(N,1);
end

% Rearranging for sinks
%NumSinksPerSource=zeros(size(Sources));
%for i=1:length(Sources)
%    NumSinksPerSource(i)=length(Sinks{i});
%end
%CumNumSinksPerSource=cumsum(NumSinksPerSource);
%TotNumSinks=CumNumSinksPerSource(end);
%CumNumSinksPerSource=[0; CumNumSinksPerSource(1:end-1)];


StateSize= E*(T-1) + (T-1)*S;
numFlowVariables = E*(T-1);


if debugflag
    fprintf('State size: %d, of which %d are flow variables \n',StateSize, numFlowVariables)
end


FindRoadLinkRtij=     @(t,i,j)   (t-1)*E + cumRoadNeighbors(i) + RoadNeighborCounter(i,j);
FindSourceRelaxti= @(t,i) (T-1)*E + (t-1)*S + i;
%FindStartRi= @(i)   T*E + i;
%FindBreaksRi= @(i)   T*E + S + i;


%% COST
if debugflag
    fprintf('Building cost: travel time...')
end
f_cost=zeros(StateSize,1);
% rebalancers' travel time
for i=1:N
    for j=RoadGraph{i}
        for t=1:(T-1)
            if i ~= j
                f_cost(FindRoadLinkRtij(t,i,j))= RebWeight*TravelTimes(i,j); %
            else
                f_cost(FindRoadLinkRtij(t,i,j))= TravelTimes(i,j);
            end
            f_cost(FindSourceRelaxti(t,i)) = SourceRelaxCost;
        end
    end
end

%% INITIALIZING CONSTRAINTS
if (debugflag)
    disp('Initializing constraints')
end

% Vehicles: N*(M+1)*C + S + TotNumSinks equality constraints, one per node and per flow and
%   per charge level plus one for each source (for the fractions), 
%   and one per class (where we choose the charge distribution at the sink of the flows)
% 2*(E+NumChargers)*(M+1)*C + 2*M*C + 2*TotNumSinks*C + C*S + C*TotNumSinks entries, two per link (including chargers) per
%   flow per charge level plus three per class, two in the conservation plus one
%   in the sum(FindPaxSinkChargeck(:,k)) = Flows(k), and one per fraction 
%   for the sum(FindChargeFractioncs(:,s)) =1

n_eq_constr = N*(T-1);
n_eq_entries = 2*E*(T-1);

if sourcerelaxflag
    n_eq_entries=n_eq_entries+S*(T-1);
end

Aeqsparse=zeros(n_eq_entries,3);
Beq=zeros(n_eq_constr,1);
Aeqrow=1;
Aeqentry=1;

% Vehicles: E inequality constraints, one per road. Each inequality
%   constraint has (M+1)*C + 1 entries one per flow per charge level, incl. reb. and one for the
%   relaxation.

if debugflag
    fprintf('Building LP program with statesize %d,%d equality constraints with %d entries', ...
        StateSize,n_eq_constr, n_eq_entries)
end


%% EQUALITY CONSTRAINTS

if (debugflag)
    disp('Building sparse equality constraints...')
end
% Conservation of rebalancers
if debugflag
    disp('Building road map for rebalancers')
    fprintf('Time step: ')
end
% we want the remaining vehicles to be evenly distributed
vdesired = max(floor((sum(Starters) + sum(sum(FlowsIn)) - sum(sum(FlowsOut))) / S), 0);
Breakers = zeros(S,1);
Breakers(1:end) = vdesired;
for t=1:T-1
    if debugflag
        fprintf(' %d/%d ',t,T)
    end
    for i=1:N
        if ~isempty(RoadGraph{i})
            for j=RoadGraph{i} %Out-flows
                if (TravelTimes(i,j)+t<=T)
                    Aeqsparse(Aeqentry,:)=[Aeqrow,FindRoadLinkRtij(t,i,j), 1];
                    Aeqentry=Aeqentry+1;
                end
            end
        end
        if ~isempty(ReverseRoadGraph{i})
            for j=ReverseRoadGraph{i} %In-flows
                if (TravelTimes(j,i) < t)
                    Aeqsparse(Aeqentry,:)=[Aeqrow,FindRoadLinkRtij(t - TravelTimes(j,i),j,i),-1];
                    Aeqentry=Aeqentry+1;
                end
            end
        end
        if sourcerelaxflag
            Aeqsparse(Aeqentry,:)=[Aeqrow,FindSourceRelaxti(t,i),-1];
            Aeqentry=Aeqentry+1;
        end
        if t==1
            Beq(Aeqrow)= Starters(i) + FlowsIn(t,i) - FlowsOut(t,i);
        elseif t == T-1
            Beq(Aeqrow)= FlowsIn(t,i) - FlowsOut(t,i) - Breakers(i);
        else
            Beq(Aeqrow)= FlowsIn(t,i) - FlowsOut(t,i);
        end
        Aeqrow=Aeqrow+1;
    end
end

if debugflag
    disp('Done! Now moving to inequalities...')
end 

%% INEQUALITY CONSTRAINTS
if debugflag
    disp('Building sparse inequality constraints...')
end

%% Make equality and inequality matrices

if Aeqrow-1~=n_eq_constr
    fprintf('ERROR: unexpected number of equality constraints (expected: %d, actual: %d)\n',n_eq_constr,Aeqrow-1)
end
if Aeqentry-1~=n_eq_entries
    fprintf('Warning: unexpected number of equality entries (expected: %d, actual: %d)\n',n_eq_entries,Aeqentry-1)
end

if (debugflag)
    disp('Building matrices from sparse representation')
end

Aeqsparse=Aeqsparse(1:Aeqentry-1,:);

Aeq=sparse(Aeqsparse(:,1),Aeqsparse(:,2),Aeqsparse(:,3), Aeqrow-1, StateSize);


%% Upper and lower bounds
if (debugflag)
    disp('Building upper and lower bounds')
end

lb=zeros(StateSize,1); %Passenger and rebalancing flows, passenger sources and sinks,...
%                       generator loads, phase angles, slack variables can't be negative
ub=Inf*ones(StateSize,1); %Why not? We enforce capacity separately





%% Call optimizer

if (debugflag)
    disp('Calling optimizer')
end
tic
options = cplexoptimset('Display', 'on');% 'Algorithm', 'interior-point');
%options = cplexoptimset('Display', 'on', 'Algorithm', 'interior-point');
%options.barrier.crossover  = -1;
%options.barrier.limits.objrange = 1e50;
[cplex_out,fval,exitflag,output]=cplexlp(f_cost,[],[],Aeq,Beq,lb,ub,[],options) ;
toc


if (debugflag)
    fprintf('Solved! fval: %f\n', fval)
    disp(output)
    %fval
    %exitflag
    %output
end

rebalanceQueue = cell(S,1);
for i = 1:S
    for j = 1:S
        if i ~= j
            for k = 1:cplex_out(FindRoadLinkRtij(1,i,j))
                rebalanceQueue{i} = [rebalanceQueue{i} j];
            end
        end
    end
end

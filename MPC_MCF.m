function [rebalanceQueue]=MPC_MCF(RoadNetwork,RebWeight,Passengers,Flags)

% Taken from AMoD-power:TIBalancedPowerFlow and edited as needed.

% Takes a road network (RoadNetwork), and a set of predicted passenger demands (Passengers) and 
% returns the first rebalancing actions of the optimal policy.

% Inputs:
% - struct RoadNetwork
% -- cell{nx1} RoadGraph: RoadGraph{i} contains the neighbors of node i in the road graph
% -- matrix(nxn) TravelTimes: TravelTimes(i,j) is the quantized travel time between nodes i and j
% -- int T: Quantized length of the planning horizon
% -- matrix(nx1) Starters: Starters(t,i) gives the number of available vehicles at station i and time t
% 
% - float RebWeight: relative cost between of idle vehicles vs traveling vehicles
%
% - struct Passengers
% -- matrix(lx3) FlowsOut: Listing of predicted flows. l is arbitrary.
% --                        FlowsOut(i,:) := [origin node, destination node, departure_time, expected demand] for the ith flow.
% -- 
% 
% - struct Flags
% -- int milpflag: if 1 solve as Mixed-Integer Linear Program (MILP)
% -- float sourcerelaxflag: if 1 add slack variables to relax the contraint of meeting demand



%% Unpack things

T=RoadNetwork.T;
RoadGraph=RoadNetwork.RoadGraph;
TravelTimes=RoadNetwork.TravelTimes;
Starters=RoadNetwork.Starters; % Starters(i) gives the number of available vehs at t=1

FlowsOut = Passengers.FlowsOut;

milpflag=Flags.milpflag;
%sourcerelaxflag=Flags.sourcerelaxflag;

%% Utilities

debugflag=1;        %Makes output verbose
DIAGNOSTIC_FLAG=0;  % Diagnoses state allocation. Useful for initial debugging.

%CongestionCost=1e3; %The cost of violating the congestion constraint
SourceRelaxCost=1e6;     % The cost of dropping a source or sink altogether
WaitTimeCost   = SourceRelaxCost / T;     %The cost per unit of time letting customers wait.

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
M=N;
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

% build finders 

StateSize= T*E + T*E + T*E + T*E;
numFlowVariables = T*E + T*E + T*E + T*E;


if debugflag
    fprintf('State size: %d, of which %d are flow variables \n',StateSize, numFlowVariables)
end


FindRoadLinkPtij= @(t,i,j) (t-1)*E +cumRoadNeighbors(i) + RoadNeighborCounter(i,j);
FindRoadLinkRtij=     @(t,i,j)  T*E + (t-1)*E + cumRoadNeighbors(i) + RoadNeighborCounter(i,j);
FindWaitingPaxtij = @(t,i,j) T*E + T*E + (t-1)*E + cumRoadNeighbors(i) + RoadNeighborCounter(i,j); % sinks are grouped by destination!
FindRealPaxtij = @(t,i,j) T*E + T*E + T*E + (t-1)*E + cumRoadNeighbors(i) + RoadNeighborCounter(i,j);


%% COST
if debugflag
    fprintf('Building cost: travel time...')
end
f_cost=zeros(StateSize,1);


for i=1:N
    for t=1:T
        for j=RoadGraph{i}
            % rebalancing costs
            if i ~= j
                f_cost(FindRoadLinkRtij(t,i,j))= RebWeight*TravelTimes(i,j); %
            else
                f_cost(FindRoadLinkRtij(t,i,j))= TravelTimes(i,j);
            end
            % waiting costs
            f_cost(FindWaitingPaxtij(t,i,j))= SourceRelaxCost;
            f_cost(FindRealPaxtij(t,i,j))= WaitTimeCost*t;
        end
    end
end

%% INITIALIZING CONSTRAINTS
if (debugflag)
    disp('Initializing constraints')
end

n_eq_constr = T*N*N + T*N + N*N;
n_eq_entries = T*N*N*3 + T*N*N*4 + N*N*T;

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
vdesired = max(floor(sum(sum(Starters)) / S), 0);
Breakers = zeros(S,1);
Breakers(1:end) = vdesired;
% there are two types of constraints, those that i) deal with pax routing, ii) deal with veh consistency

for t=1:T
    if debugflag
        fprintf(' %d/%d ',t,T)
    end
    
    for i=1:N
        % i) pax constraints
        if ~isempty(RoadGraph{i})
            for j=RoadGraph{i} %Out-flows, pax embark to their destination
                %if (TravelTimes(i,j)+t<=T)
                    Aeqsparse(Aeqentry,:)=[Aeqrow,FindRoadLinkPtij(t,i,j), 1];
                    Aeqentry=Aeqentry+1;
                %end
                % waiting, we get to postpone the departure
                Aeqsparse(Aeqentry,:)=[Aeqrow,FindWaitingPaxtij(t,i,j), 1];
                Aeqentry=Aeqentry+1;
                % customers that waited the last timestep
                %if t > 1
                %    Aeqsparse(Aeqentry,:)=[Aeqrow,FindWaitingPaxtij(t-1,i,j), -1];
                %    Aeqentry=Aeqentry+1;
                %end
                % initial customers that get to wait
                Aeqsparse(Aeqentry,:)=[Aeqrow,FindRealPaxtij(t,i,j), -1];
                Aeqentry=Aeqentry+1;
                % all of this must equal the predicted demand for (t,i,j)
                % and we want this to equal the total amount of demand at that time and node
                if t > 1
                    Beq(Aeqrow)= FlowsOut{t}(i,j);
                else
                    Beq(Aeqrow) = 0;
                end
                Aeqrow=Aeqrow+1;
            end
        end

        % ii) vehicle constraints
        if ~isempty(RoadGraph{i})
            for j=RoadGraph{i} %Out-flows
                %if (TravelTimes(i,j)+t<=T)
                    % rebalancing outlfows
                    Aeqsparse(Aeqentry,:)=[Aeqrow,FindRoadLinkRtij(t,i,j), 1];
                    Aeqentry=Aeqentry+1;
                    % pax outflows
                    Aeqsparse(Aeqentry,:)=[Aeqrow,FindRoadLinkPtij(t,i,j), 1];
                    Aeqentry=Aeqentry+1;
                %end
            end
        end
        if ~isempty(ReverseRoadGraph{i})
            for j=ReverseRoadGraph{i} %In-flows
                if (TravelTimes(j,i) < t)
                    % rebalancing inflows
                    Aeqsparse(Aeqentry,:)=[Aeqrow,FindRoadLinkRtij(t - TravelTimes(j,i),j,i),-1];
                    Aeqentry=Aeqentry+1;
                    % pax inflows
                    Aeqsparse(Aeqentry,:)=[Aeqrow,FindRoadLinkPtij(t - TravelTimes(j,i),j,i), -1];
                    Aeqentry=Aeqentry+1;
                end
            end
        end
        if t == T
            Beq(Aeqrow)= Starters(t,i);% - Breakers(i);
        else
            Beq(Aeqrow)= Starters(t,i);
        end
        Aeqrow=Aeqrow+1;
    end
end

for i=1:N
    for j=1:N
        for t=1:T
            Aeqsparse(Aeqentry,:)=[Aeqrow,FindRealPaxtij(t,i,j), 1];
            Aeqentry=Aeqentry+1;
        end
        Beq(Aeqrow)= FlowsOut{1}(i,j);
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

lb=zeros(StateSize,1); %everything is non-negative
ub=Inf*ones(StateSize,1); %no constraints





%% Call optimizer

if (debugflag)
    disp('Calling optimizer')
end
%tic
%options = cplexoptimset('Display', 'on');% 'Algorithm', 'interior-point');
%%options = cplexoptimset('Display', 'on', 'Algorithm', 'interior-point');
%%options.barrier.crossover  = -1;
%%options.barrier.limits.objrange = 1e50;
%[cplex_out,fval,exitflag,output]=cplexlp(f_cost,[],[],Aeq,Beq,lb,ub,[],options) ;
%toc


%% Variable type

if (debugflag)
    disp('Building constraint type')
end

% Building block
ConstrType=char(zeros(1,StateSize));

% Continuous-probability version
if (~milpflag)
    ConstrType(1:end)='C';
else
    % Actual MILP formulation
    ConstrType(1:end)='I';
end


if (milpflag)
    fprintf('Solving as MILP.')
    MyOptions=cplexoptimset('cplex');;
    %MyOptions.parallel=1;
    %MyOptions.threads=8;
    %MyOptions.mip.tolerances.mipgap=0.01;
    %MyOptions.Display = 'iter';
    if (debugflag)
        tic
    end
    [cplex_out,fval,exitflag,output]=cplexmilp(f_cost,[],[],Aeq,Beq,[],[],[],lb,ub,ConstrType,[],MyOptions);
    
    if (debugflag)
        toc
    end
else
    if (debugflag)
        tic
    end
    [cplex_out,fval,exitflag,output]=cplexlp(f_cost,[],[],Aeq,Beq,lb,ub) ;
    if (debugflag)
        toc
    end
end


if (debugflag)
    fprintf('Solved! fval: %f\n', fval)
    disp(output)
    %fval
    %exitflag
    %output
end

delivered = 0;
dropped = 0;
waiting = 0;
total_wait = 0;
for t=1:T
    for i=1:N
        for j=1:N
            delivered = delivered + cplex_out(FindRoadLinkPtij(t,i,j));
            dropped = dropped + cplex_out(FindWaitingPaxtij(t,i,j));
            if t > 1
                w = cplex_out(FindRealPaxtij(t,i,j));
                waiting = waiting + w;
                if w > 0
                    total_wait = total_wait + w*(t-1);
                end
            end
        end
    end
end
fprintf('Total delivered: %f\n', delivered)
fprintf('Total dropped: %f\n', dropped)
fprintf('Total waiting: %f\n', waiting)
fprintf('Average waiting: %f\n', total_wait / waiting)

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

% Tests MPC with waiting
RoadGraph = {};
for i=1:5
    RoadGraph{i} = [1 2 3 4 5];
end

TravelTimes = [ 1 1 11 21 22 ; 1 1 10 20 21; 11 10 1 10 11; 21 20 10 1 1; 22 21 11 1 1];


v = 1;
N = length(RoadGraph);

T = 50;
RoadNetwork.T = T;
RoadNetwork.RoadGraph = RoadGraph;
RoadNetwork.TravelTimes = TravelTimes;

Flags.milpflag = 0;

RebWeight = 5.0;

Starters = zeros(T,length(RoadGraph));
Starters(1,3) = 1

RoadNetwork.Starters = Starters;

% random demand
FlowsOut = cell(1,T);
counter = 1;
tot_pax = 0;
for t=1:T
    FlowsOut{t} = sparse(N,N);
end
FlowsOut{1}(2,1) = 1;
FlowsOut{1}(4,5) = 1;

Passengers.FlowsOut = FlowsOut;

[rebalanceQueue, cplex_out]=MPC_MCF(RoadNetwork,RebWeight,Passengers,Flags);


unique(abs(cplex_out -round(cplex_out)))

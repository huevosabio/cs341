% Tests MPC with waiting
BuildCityMap;
load('didi.mat');
v = 5000;
N = length(RoadGraph);

T = 50;
RoadNetwork.T = T;
RoadNetwork.RoadGraph = RoadGraph;
RoadNetwork.TravelTimes = TravelTimes;

Flags.milpflag = 0;

RebWeight = 5.0;

Starters = zeros(T,length(RoadGraph));

% placing the vehicles randomly
Tinit = 25;
R = rand(v,2);
R(:,1) = ceil(R(:,1) * N);
R(:,2) = ceil(R(:,2) .^2 * Tinit);

for t=1:Tinit
    for i=1:N
        mask = R(:,1) == i & R(:,2) == t;
        Starters(t,i) = sum(mask);
    end
end
RoadNetwork.Starters = Starters;

% random demand
FlowsOut = cell(1,T);
counter = 1;
tot_pax = 0;
for t=1:T
    FlowsOut{t} = sparse(N,N);
    for i=1:N
        for j=1:N
            if t == 1
                FlowsOut{t}(i,j) =  poissrnd(0.5);
                tot_pax = tot_pax + FlowsOut{t}(i,j);
            end
        end
    end
end

Passengers.FlowsOut = FlowsOut;

[rebalanceQueue]=MPC_MCF(RoadNetwork,RebWeight,Passengers,Flags);

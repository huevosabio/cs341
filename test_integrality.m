ntrials = 100;

% Tests MPC with waiting
BuildCityMap;
load('didi.mat');
v = 5000;
N = length(RoadGraph);

Tinit = 25;
R = rand(v,2);
R(:,1) = ceil(R(:,1) * N);
R(:,2) = ceil(R(:,2) .^2 * Tinit);

T = 50;
RoadNetwork.T = T;
RoadNetwork.RoadGraph = RoadGraph;
RoadNetwork.TravelTimes = TravelTimes;

Flags.milpflag = 0;

RebWeight = 5.0;

Starters = zeros(T,length(RoadGraph));

nints = 0;

results = zeros(ntrials,8);

for k=1:ntrials
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
	            if t <= Tinit
	                FlowsOut{t}(i,j) =  poissrnd(0.2);
	                tot_pax = tot_pax + FlowsOut{t}(i,j);
	            end
	        end
	    end
	end
	Passengers.FlowsOut = FlowsOut;
	Flags.milpflag = 0;
	[lrebalanceQueue, lintegral, output, lfval, lrfval, lfval2, lrfval2, lrconstr]=MPC_MCF(RoadNetwork,RebWeight,Passengers,Flags);
	lptime = output.time;
	Flags.milpflag = 1;
	% 
	[mrebalanceQueue, mintegral, output, fval, rfval, mfval2, mrfval2, rconstr]=MPC_MCF(RoadNetwork,RebWeight,Passengers,Flags);
	milptime = output.time;
	nints = nints + lintegral;
	ig = fval / lfval;
	rig = fval / lrfval;
	ig2 = mfval2 / lfval2;
	rig2 = mfval2 / lrfval2;
	results(k,:) = [lptime, milptime, lintegral, ig, rig, ig2, rig2, lrconstr];
	fprintf('number of non-integrals: %d / %d \n', nints, k)
end
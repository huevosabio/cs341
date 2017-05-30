% runs the didi problem
clear all; close all; clc;
load('didi.mat')

RoadNetwork.T = T;
RoadNetwork.RoadGraph = RoadGraph;
%RoadNetwork.RoadCap = RoadCap;
RoadNetwork.TravelTimes = TravelTimes;

Passengers.FlowsOut = FlowsOut;
Passengers.FlowsIn = FlowsIn;

Flags.milpflag = 0;
Flags.congrelaxflag = 0;
Flags.sourcerelaxflag = 0;
Flags.cachedAeqflag = 0;

RebWeight = 5.0;
[cplex_out, Data, FindFunctions]=SolveRebalancing(RoadNetwork,RebWeight,Passengers,Flags);

RebIn = Data.RebIn;
RebOut = Data.RebOut;
StayIn = Data.StayIn;
StayOut = Data.StayOut;
nvehs = Data.nvehs;

save('Data', 'RebIn','RebOut', 'StayIn', 'StayOut', 'nvehs')    

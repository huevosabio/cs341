%LinkTime=csvread('travel_time.csv',4,1);
%LinkTime=LinkTime*60; %express in seconds

global RoadGraph NodesLocation

NodesLocation=csvread('inferred_locations.csv',1,1);
N=length(NodesLocation);

NodesLocation=NodesLocation*60;

StationLocation=NodesLocation;
StationNodeID=[1:1:N]';


RoadGraph=cell(N,1);
LinkFreeFlow=zeros(N,N);
LinkLength=LinkFreeFlow;
RoadCap=LinkFreeFlow;
MAXCAP=1e6;

for i=1:N
    for j=1:N
        RoadGraph{i}=[RoadGraph{i},j];
        LinkLength(i,j)=norm(NodesLocation(j,:) - NodesLocation(i,:));
        LinkTime(i,j)=LinkLength(i,j) * 60; % time in seconds, lengths are in minutes! - RI
        LinkFreeFlow(i,j)=LinkLength(i,j)/LinkTime(i,j);
        RoadCap(i,j)=MAXCAP;
    end
end
LinkTime(LinkTime==0)=1e9;
LinkLength(LinkLength==0)=1e9;

LinkNumVehicles = sparse(N,N);
LinkCapacityLeft = sparse(N,N);
LinkTime=csvread('travel_time.csv',4,1);

global RoadGraph NodesLocation

NodesLocation=csvread('inferred_locations.csv',1,1);
N=length(NodesLocation);



StationLocation=NodesLocation;
StationNodeID=[1:1:N]';


RoadGraph=cell(N,1);
LinkFreeFlow=zeros(N,N);
LinkLength=LinkFreeFlow;
RoadCap=LinkFreeFlow;
MAXCAP=1e6;

for i=1:N
    for j=1:N
        if LinkTime(i,j)>0
            RoadGraph{i}=[RoadGraph{i},j];
            LinkLength(i,j)=norm(NodesLocation(j,:) - NodesLocation(i,:));
            LinkFreeFlow(i,j)=LinkLength(i,j)/LinkTime(i,j);
            RoadCap=MAXCAP;
        end
    end
end

LinkNumVehicles = sparse(N,N);
LinkCapacityLeft = sparse(N,N);
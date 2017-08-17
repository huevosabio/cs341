function [Route] = findRoute(start,goal, LinkTime)
global RoadGraph NodesLocation;
ClosedSet = [];
OpenSet = [start];

N=size(RoadGraph,1);

% Heuristic, specific to this square grid implementation

persistent heuristic_est
if isempty(heuristic_est)
    heuristic_est=zeros(N,N);
    for ii=1:N
        for jj=1:N
            heuristic_est(ii,jj) = norm(NodesLocation(jj,:)-NodesLocation(ii,:), 2);
        end
    end
end
% /heuristic

CameFrom=-1*ones(N);
g_score = Inf*ones(N,1);
f_score = Inf*ones(N,1);
g_score(start)=0;

f_score(start)=g_score(start)+heuristic_est(start,goal);

while ~isempty(OpenSet)
    [~,currInOpenSet]=min(f_score(OpenSet));
    currNode = OpenSet(currInOpenSet);
    if currNode == goal
        
        Route = goal;
        while (Route(1)~=start)
            Route=[CameFrom(Route(1)),Route];
        end
        return
    end
    OpenSet = OpenSet(OpenSet~=currNode);       %Remove current node
    ClosedSet = [ClosedSet, currNode];          %Add to closed set
    for neigh = RoadGraph{currNode}
        if ~sum(ClosedSet == neigh)             %If neighbor is not in closed set
            tentative_g = g_score(currNode) + LinkTime(currNode, neigh); 
            if (~sum(neigh == OpenSet))         %If the neighbor is not already in OpenSet, add it
                OpenSet = [OpenSet, neigh];     
            elseif tentative_g>=g_score(neigh)
                continue
            end
            CameFrom(neigh) = currNode;
            g_score(neigh) = tentative_g;
            f_score(neigh) = g_score(neigh)+heuristic_est(neigh,goal);
            
        end
    end
end


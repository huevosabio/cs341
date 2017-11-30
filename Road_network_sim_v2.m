% Ported from AMoD_congestion on June 1, 2017

%road map specifications:
% - list of node positions
% - list of neighbors for each node
% - link capacities, where each link is defined by the tuple (u,v). The
% link capacities will be implemented as a sparse matrix
% - link freeflow speed

% Congestion model
% travel time function using the BPR formula, gives the travel time as a
% function of how many vehicles are on the link. 
% Option 1: The speed of each link at 
% each time step is determined by the number of vehicles in that link at
% that time. ---> *use this one*
% **Option 2: The travel time of each
% vehicle is determined upon its arrival to a link. The vehicle's speed is
% then found from its travel time. 
% 1. group all vehicles entering a link.
% 2. calculate speeds for all vehicles that just entered the link


% Simulation parameters:
% 1. for each link (u,v)
% - unit vector for the direction of each link
% - number of vehicles in the link
% - current speed for the link (calculated every time step based on the
% utilization/number of vehicles)

% 2. for each vehicle
% - path from origin to destination (path calculated with A*). Path should
% be a list of nodes
% - current position vector
% - current stage of trip (index of which node it most recently passed).
% This also provides the current link the vehicle is on

% 3. 

clear all;
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio126\cplex\matlab\x64_win64');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeStep = 10;      % number of time steps per minute
dt = 60/timeStep;   % length of each time step
Tmax = 60*24*timeStep;    % simulation time
Thor = 5*timeStep;      % rebalancing time horizon (for real-time algorithm) for vehicles
Tthresh = 0;     % at what time to start gathering data
v = 5000;
rebPaths = cell(0);
ccTmp = 1;
cc = 1;

PLOTFLAG = 0;
PLOTREBFLAG = 1;
NAIVEFLAG = 1;
rebCount = 1;
MPCFLAG = 1;

SIMNAME = 'MPC-EMPTY';

if MPCFLAG
    % HACKY SHIT: LOADS THE REAL TRIP COUNT BY 5MIN
    load('didi.mat')
    load('ignored_assets/tod_predictions_empty.mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BuildCityMap; 
% N = number of nodes
% RoadGraph = list of neighbors
% NodesLocation = x and y location of each node
% NodesLatLon = lat and lon of each node
% 


if MPCFLAG

    % MPC variables
    predictionStep = 5*timeStep; % number of steps a prediction covers
    horizon = 50; % longest travel time is 44 steps...

    RoadNetwork.T = horizon;
    RoadNetwork.RoadGraph = RoadGraph;
    RoadNetwork.TravelTimes = TravelTimes;

    Flags.milpflag = 1;
    Flags.congrelaxflag = 0;
    Flags.sourcerelaxflag = 1;
    Flags.cachedAeqflag = 0;

    RebWeight = 5.0; % how much more expensive it is to rebalance than to remain idle

end

LoadStateDefinitions;

numStations = length(StationNodeID);

for i = 1:numStations
    station(i) = struct('id',i,'carIdle',[], 'custId',[], 'carOnRoad',[],...
        'custUnassigned',[],'waitTimes',[], 'arriveHour', [], 'node_id', StationNodeID(i));
end
% for naive rebalancing
if NAIVEFLAG
    rebalanceQueue = cell(numStations,1);
    % define Tij
    %Tij = zeros(numStations, numStations);
    %for i = 1:numStations
    %    for j = 1:numStations
    %        Tij(i,j) = norm(StationLocation(j,:) - StationLocation(i,:), 1);
    %    end
    %end
    Tij=LinkTime;
end

% data collection
avgWaitTimes = zeros(Tmax,1);
cumNumCustomers = zeros(Tmax,1);
numVehiclesBusy = zeros(Tmax, 1);
numVehiclesRebalancing = zeros(Tmax, 1);
numVehiclesDrivingToPickup = numVehiclesRebalancing;
numVehiclesDrivingToDest = numVehiclesRebalancing;
numVehiclesDrivingToStation = numVehiclesRebalancing;
numVehiclesOnSelfLoop =  numVehiclesRebalancing;
numVehiclesIdle = numVehiclesRebalancing;
numVehiclesNotRebalancing = zeros(Tmax, 1);
%numCarsOnLink = cell(Tmax);
carsIdlePrint=zeros(Tmax,numStations);
carsOnRoadPrint=zeros(Tmax,numStations);
custUnassigned=zeros(Tmax,numStations);
custWaiting=zeros(Tmax,numStations);
%
numRebTasks = zeros(Tmax,numStations);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load demand data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Loading demand data...')
filename = 'ignored_assets/MATLAB_orders.csv';
MData = csvread(filename,1,1);
arrivalTimeOffset = (0*3600 + -5*1*60 + 0)/60*timeStep;
arrivalTimes = (MData(:,3)*3600 + MData(:,4)*60 + MData(:,5))/60*timeStep - arrivalTimeOffset;
fprintf('loaded!\n')

% but just start at the arrivalTimeOffset
while arrivalTimes(ccTmp) < 0
    ccTmp = ccTmp + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% declare data structures for cars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% struct for each car
stationCounter = 1;
for i = 1:v
    car(i) = struct('id',i,'passId', 0, 'dstation',stationCounter,'ostation',stationCounter,...
        'dpos', [], 'state', IDLE, 'pos', NodesLocation(station(stationCounter).node_id, :),...
        'direction',[], 'path', [StationNodeID(stationCounter)], 'time_left', 0, 'speedfactor', 1);
    % update station data for this car
    station(stationCounter).carIdle = [station(stationCounter).carIdle, i];
    if stationCounter == numStations
        stationCounter = 1;
    else
        stationCounter = stationCounter + 1;
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin Simulation: main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:Tmax-1
    % update LinkTime variable
    % go through all the links (i,j)
    
    %No congestion
    %for i = 1:N
    %    for j = RoadGraph{i}
    %        [tmpLinkTime, tmpLinkSpeed] = getSpeed(LinkNumVehicles(i,j), RoadCap(i,j), LinkFreeFlow(i,j), LinkLength(i,j));
    %        LinkTime(i,j) = tmpLinkTime;
    %        LinkSpeed(i,j) = tmpLinkSpeed;
    %    end
    %end
    LinkSpeed=LinkFreeFlow;
    
    for i=1:numStations
        carsIdlePrint(t,i)=length(station(i).carIdle);
        carsOnRoadPrint(t,i)=length(station(i).carOnRoad);
        custUnassigned(t,i)=length(station(i).custUnassigned);
        custWaiting(t,i)=length(station(i).custUnassigned) + length(station(i).custId);
    end
    carsIdlePrint(t,:)+carsOnRoadPrint(t,:);

    fprintf('Time: %d, Unassigned Customers: %d, Waiting Customers: %d, Total Customers: %d \n', t, sum(custUnassigned(t,:)), sum(custWaiting(t,:)), cc-1)
    % vehicle state transitions
    
    for i = 1:v
        if car(i).state ~= IDLE && car(i).state ~= SELF_LOOP && length(car(i).path) == 2 % if it's on the last leg of its trip
            distToDest = norm(car(i).dpos - car(i).pos,2);
            if distToDest < car(i).speedfactor * LinkSpeed(car(i).path(1), car(i).path(2))*dt
                % remove car from link
                LinkNumVehicles(car(i).path(1), car(i).path(2)) = LinkNumVehicles(car(i).path(1), car(i).path(2)) - 1;
                % shrink path to just final destination
                car(i).path = car(i).path(2);
                
                % Blurb: if the car is on a self-loop and got to the
                % destination, release it as if it was driving_to_dest. If
                % the car is driving to pick up and gets a customer and the
                % customer is on a self-loop, set travel_time and put in
                % status self_loop. Actually, car never drive to a
                % customer, since each node is a station: scratch that, go
                % easy.
                
                if car(i).state == DRIVING_TO_PICKUP
                    % pick up customer
                    if customer(car(i).passId).onode == customer(car(i).passId).dnode
                        car(i).state = SELF_LOOP;
                        car(i).dpos = customer(car(i).passId).dpos;
                        car(i).dstation = customer(car(i).passId).dstation;
                        car(i).ostation = customer(car(i).passId).ostation;
                        car(i).time_left = customer(car(i).passId).traveltime * 60;

                        customer(car(i).passId).pickedup = 1;
                        % remove customer from the origin station
                        tmpindex = find(station(customer(car(i).passId).ostation).custId == car(i).passId);
                        station(customer(car(i).passId).ostation).custId = [station(customer(car(i).passId).ostation).custId(1:tmpindex-1), station(customer(car(i).passId).ostation).custId(tmpindex+1:end)];
                    else
                        car(i).state = DRIVING_TO_DEST;
                        car(i).pos = car(i).dpos;
                        if isnan(car(i).pos)
                            'Error: car position is NaN'
                        end
                        car(i).dpos = customer(car(i).passId).dpos;
                        customer(car(i).passId).pickedup = 1;
                        % remove customer from the origin station
                        tmpindex = find(station(customer(car(i).passId).ostation).custId == car(i).passId);
                        station(customer(car(i).passId).ostation).custId = [station(customer(car(i).passId).ostation).custId(1:tmpindex-1), station(customer(car(i).passId).ostation).custId(tmpindex+1:end)];
                        %
                        car(i).path = findRoute(car(i).path(1), customer(car(i).passId).dnode, LinkTime);
                        car(i).speedfactor =  LinkTime(customer(car(i).passId).onode,customer(car(i).passId).dnode) / (customer(car(i).passId).traveltime * 60);
                        
                        LinkNumVehicles(car(i).path(1), car(i).path(2)) = LinkNumVehicles(car(i).path(1), car(i).path(2)) + 1;
                        car(i).direction = (NodesLocation(car(i).path(2),:) - NodesLocation(car(i).path(1),:)); car(i).direction=car(i).direction/norm(car(i).direction);
                    end
                elseif car(i).state == DRIVING_TO_DEST
                    % drop off customer
                    car(i).state = DRIVING_TO_STATION;
                    car(i).pos = customer(car(i).passId).dpos;
                    if isnan(car(i).pos)
                        'Error: car position is NaN'
                    end
                    car(i).dpos = StationLocation(customer(car(i).passId).dstation, :);
                    % route car to station
                    car(i).path = findRoute(car(i).path(1), station(customer(car(i).passId).dstation).node_id, LinkTime);
                    % car is now free to receive assignments from station
                    
                    customer(car(i).passId).delivered = 1;
                    car(i).passId = 0;

                    % return to normal speed
                    car(i).speedfactor = 1;
                    
                    % if the destination was already a station
                    if length(car(i).path) == 1
                        car(i).state = IDLE;
                        car(i).ostation = car(i).dstation;
                        station(car(i).ostation).carIdle = [station(car(i).ostation).carIdle, car(i).id];
                    else
                        station(car(i).dstation).carOnRoad = [station(car(i).dstation).carOnRoad, car(i).id];
                        LinkNumVehicles(car(i).path(1), car(i).path(2)) = LinkNumVehicles(car(i).path(1), car(i).path(2)) + 1;
                        car(i).direction = (NodesLocation(car(i).path(2),:) - NodesLocation(car(i).path(1),:)); car(i).direction=car(i).direction/norm(car(i).direction);
                    end
                    
                elseif car(i).state == DRIVING_TO_STATION
                    % arrived back to station
                    car(i).state = IDLE;
                    car(i).pos = car(i).dpos;
                    if isnan(car(i).pos)
                        'Error: car position is NaN'
                    end
                    car(i).ostation = car(i).dstation;
                    % add car to idle car list
                    station(car(i).ostation).carIdle = [station(car(i).ostation).carIdle, car(i).id];
                    tmpindex = find(station(car(i).ostation).carOnRoad == car(i).id);
                    % remove car from caronroad list
                    station(car(i).ostation).carOnRoad = [station(car(i).ostation).carOnRoad(1:tmpindex-1), station(car(i).ostation).carOnRoad(tmpindex+1:end)];
                    
                elseif car(i).state == REBALANCING
                    % finished rebalancing to station
                    car(i).state = IDLE;
                    car(i).pos = car(i).dpos;
                    if isnan(car(i).pos)
                        'Error: car position is NaN'
                    end
                    car(i).ostation = car(i).dstation;
                    % add car to idle car list
                    station(car(i).ostation).carIdle = [station(car(i).ostation).carIdle, car(i).id];
                    tmpindex = find(station(car(i).ostation).carOnRoad == car(i).id);
                    % remove car from caronroad list
                    if ~isempty(tmpindex)
                        station(car(i).ostation).carOnRoad = [station(car(i).ostation).carOnRoad(1:tmpindex-1), station(car(i).ostation).carOnRoad(tmpindex+1:end)];
                    end
                end
            end
        end
        if car(i).state == SELF_LOOP && car(i).time_left < dt
            car(i).state = DRIVING_TO_STATION;
            car(i).pos = customer(car(i).passId).dpos;
            if isnan(car(i).pos)
                'Error: car position is NaN'
            end
            car(i).dpos = StationLocation(customer(car(i).passId).dstation, :);
            % route car to station
            car(i).path = findRoute(car(i).path(1), station(customer(car(i).passId).dstation).node_id, LinkTime);
            % car is now free to receive assignments from station
            
            customer(car(i).passId).delivered = 1;
            car(i).passId = 0;
            
            % return to normal speed
            car(i).speedfactor = 1;
            
            % if the destination was already a station
            if length(car(i).path) == 1
                car(i).state = IDLE;
                car(i).ostation = car(i).dstation;
                station(car(i).ostation).carIdle = [station(car(i).ostation).carIdle, car(i).id];
            else
                station(car(i).dstation).carOnRoad = [station(car(i).dstation).carOnRoad, car(i).id];
                LinkNumVehicles(car(i).path(1), car(i).path(2)) = LinkNumVehicles(car(i).path(1), car(i).path(2)) + 1;
                car(i).direction = (NodesLocation(car(i).path(2),:) - NodesLocation(car(i).path(1),:)); car(i).direction=car(i).direction/norm(car(i).direction);
            end
            
        end
    end
    
    % new customer arrivals
    if ccTmp <= max(size(arrivalTimes))
        while arrivalTimes(ccTmp) < t
            %customer arrival and destination locations
            tmpCust = [MData(ccTmp,6:7); MData(ccTmp,8:9)]*60;
            % find the nearest nodes
            tmpNodes = dsearchn(NodesLocation, tmpCust);
            %if tmpNodes(1) ~= tmpNodes(2)
            % find the stations
            tmpStations = dsearchn(StationLocation, tmpCust);
            % make customer structure
            customer(cc) = struct('opos',NodesLocation(tmpNodes(1),:), 'dpos', NodesLocation(tmpNodes(2),:),...
                'onode', tmpNodes(1), 'dnode', tmpNodes(2), 'ostation',tmpStations(1), 'dstation', tmpStations(2), ...
                'waitTime',0,'serviceTime',0,'pickedup',0,'delivered',0, 'traveltime', MData(ccTmp,10));
            % add this customer to the station
            station(customer(cc).ostation).custId = [station(customer(cc).ostation).custId  cc];
            station(customer(cc).ostation).custUnassigned = [station(customer(cc).ostation).custUnassigned  cc];
            cc = cc + 1;
            ccTmp = ccTmp + 1;
            %else
            %    ccTmp = ccTmp + 1;
            %end
            if ccTmp > max(size(arrivalTimes))
                break
            end
        end
    end
    % assign vehicles to new customers
    for i = 1:numStations
        % as long as there are cars left or customers waiting
        while ~isempty(station(i).custUnassigned) && (~isempty(station(i).carIdle) || ~isempty(station(i).carOnRoad))
            custInd = station(i).custUnassigned(1);
            % get distance from customer 
            
            distToStation = norm(customer(custInd).opos - StationLocation(i,:));
            
            % get distance to nearest car on road
            if ~isempty(station(i).carIdle)
                assignedCarID = station(i).carIdle(1);
                shortestDistance = distToStation;
            else
                shortestDistance = inf;
            end
            carOnRoadIndex = 0;
            if ~isempty(station(i).carOnRoad)
                for j = 1:length(station(i).carOnRoad)
                    distToCar = norm(customer(custInd).opos - car(station(i).carOnRoad(j)).pos);
                    if distToCar < shortestDistance
                        shortestDistance = distToCar;
                        assignedCarID = station(i).carOnRoad(j);
                        carOnRoadIndex = j;
                    end
                end
            end
            if shortestDistance == inf
                %'ERROR: inf shortest distance (could not find a vehicle for the customer)'
            end
            % assign pax to car
            car(assignedCarID).passId = custInd;
            
            tmpPath = car(assignedCarID).path;
            car(assignedCarID).path = findRoute(tmpPath(1), customer(custInd).onode, LinkTime); 
            alreadyAtPickup = 0;
            if length(tmpPath) == 1 && length(car(assignedCarID).path) == 1
                % TODO: APPLIES TO SELFLOOP
                % vehicle is not initially in motion (at a station)
                % go straight to routing delivery
                alreadyAtPickup = 1;
                if customer(custInd).onode==customer(custInd).dnode
                    car(assignedCarID).state = SELF_LOOP;
                    car(assignedCarID).dpos = customer(custInd).dpos;
                    car(assignedCarID).time_left = customer(custInd).traveltime * 60;
                else
                    car(assignedCarID).state = DRIVING_TO_DEST;
                    car(assignedCarID).dpos = customer(custInd).dpos;
                    % route the new path using the customer speed
                    car(assignedCarID).path = findRoute(tmpPath(1), customer(custInd).dnode, LinkTime);
                    car(assignedCarID).speedfactor =  LinkTime(customer(custInd).onode,customer(custInd).dnode) / (customer(custInd).traveltime * 60);
                end
                customer(custInd).pickedup = 1;
                % remove customer from the origin station
                tmpindex = find(station(customer(custInd).ostation).custId == car(assignedCarID).passId);
                station(customer(car(assignedCarID).passId).ostation).custId = [station(customer(car(assignedCarID).passId).ostation).custId(1:tmpindex-1), station(customer(car(assignedCarID).passId).ostation).custId(tmpindex+1:end)];
            else
                % TODO: APPLIES TO SELFLOOP
                % vehicle is in motion.
                if length(car(assignedCarID).path) == 1
                    car(assignedCarID).path = [tmpPath(2) tmpPath(1)];
                elseif length(tmpPath) > 1
                    if car(assignedCarID).path(2) ~= tmpPath(2)
                        % go a different way - make U-turn
                        car(assignedCarID).path = [tmpPath(2) car(assignedCarID).path];
                    end
                end
            end
            % TODO: APPLIES TO SELFLOOP
            if length(car(assignedCarID).path) > 1
                LinkNumVehicles(car(assignedCarID).path(1), car(assignedCarID).path(2)) = LinkNumVehicles(car(assignedCarID).path(1), car(assignedCarID).path(2))+1;
            end
            if length(tmpPath) > 1
                LinkNumVehicles(tmpPath(1), tmpPath(2)) = LinkNumVehicles(tmpPath(1), tmpPath(2))-1;
            end
            % TODO: APPLIES TO SELFLOOP
            % assign the car to the customer
            car(assignedCarID).dstation = customer(custInd).dstation;
            car(assignedCarID).ostation = i;
            if ~alreadyAtPickup
                car(assignedCarID).state = DRIVING_TO_PICKUP;
                car(assignedCarID).dpos = customer(custInd).opos;
            end
            % direction vector
            if length(car(assignedCarID).path) > 1
                tmpNode2 = car(assignedCarID).path(2);
                tmpNode1 = car(assignedCarID).path(1);
                car(assignedCarID).direction = (NodesLocation(tmpNode2,:) - NodesLocation(tmpNode1,:)); car(assignedCarID).direction=car(assignedCarID).direction/norm(car(assignedCarID).direction);
            end
            % remove assigned car from station
            if ~isempty(station(i).carIdle) && assignedCarID == station(i).carIdle(1)
                station(i).carIdle = station(i).carIdle(2:end);
            else
                station(i).carOnRoad = [station(i).carOnRoad(1:carOnRoadIndex-1), station(i).carOnRoad(carOnRoadIndex+1:end)];
            end
            % remove unassigned customer from station
            station(i).custUnassigned = station(i).custUnassigned(2:end);
        end
    end
    
    % update global road network capacity matrix
    
    % rebalancing
    if rebCount >= Thor
        rebCount = 1;
        
        vown = zeros(numStations,1);
        vexcess = zeros(numStations,1);
        vdesired = zeros(numStations,1);
        totalCustomers = 0;
        
        if NAIVEFLAG
            % clear rebalanceQueue
            for i = 1:numStations
                rebalanceQueue{i} = [];
            end
        end
        
        for i = 1:numStations
            vown(i) = length(station(i).carIdle) + length(station(i).carOnRoad);
        end
        for i = 1:v
            % if the car is in the process of pickup or dropoff
            if car(i).state == DRIVING_TO_DEST  || car(i).state == DRIVING_TO_PICKUP || car(i).state == SELF_LOOP
                vown(car(i).dstation) = vown(car(i).dstation) + 1;
            end
                
        end
        
        % find excess vehicles of each station
        for i = 1:numStations
            vexcess(i) = vown(i) - length(station(i).custUnassigned);
            % find total customers
            totalCustomers = totalCustomers + length(station(i).custUnassigned);
        end

        % vehicles desired for each station
        vdesired = floor((sum(vown) - totalCustomers)/numStations)*ones(numStations,1);
        
        if NAIVEFLAG
            if ~MPCFLAG
                numStations2=numStations;
                % car optimization
                cvx_begin
                    variable numij2(numStations2,numStations2)
                    minimize (sum(sum(Tij.*numij2)));
                    subject to
                    vexcess + sum((numij2' - numij2)')' >= vdesired;
                    % sum(numij')' <= rown;
                    % trace(numij) == 0;
                    numij2 >= 0;
                cvx_end
                % make sure numij is integer
                numij2 = round(numij2);            
                % add rebalancing vehicles to queues
                for i = 1:numStations
                    for j = 1:numStations
                        for k = 1:numij2(i,j)
                            rebalanceQueue{i} = [rebalanceQueue{i} j];
                            numRebTasks(t,i) = numRebTasks(t,i) + 1;
                        end
                    end
                end
            else
                % MPC stuff

                % load passenger predictions into FlowsOut + FlowsIn
                % predictions are in 5 min intervals
                currentTime = t / predictionStep;
                Passengers.FlowsOut = predictDemand(currentTime, horizon, numStations, predictor);
                RoadNetwork.Starters = zeros(horizon,numStations);
                % load current trips/rebalancing into Starters
                for i = 1:v
                    % if the car is in the process of pickup or dropoff
                    if car(i).state == DRIVING_TO_DEST  || car(i).state == DRIVING_TO_STATION || car(i).state == REBALANCING
                        distToDest = norm(car(i).dpos - car(i).pos,2);
                        eta = distToDest / car(i).speedfactor;
                        tin = ceil(eta / (predictionStep * dt));
                        if tin < horizon && tin > 0
                            RoadNetwork.Starters(tin,car(i).dstation) = RoadNetwork.Starters(tin,car(i).dstation) + 1;
                        elseif tin < horizon
                            RoadNetwork.Starters(1,car(i).dstation) = RoadNetwork.Starters(1,car(i).dstation) + 1;
                        end
                    elseif car(i).state == SELF_LOOP
                        eta = car(i).time_left;
                        tin = ceil(eta / (predictionStep * dt));
                        if tin < horizon && tin > 0
                            RoadNetwork.Starters(tin,car(i).dstation) = RoadNetwork.Starters(tin,car(i).dstation) + 1;
                        elseif tin < horizon
                            RoadNetwork.Starters(1,car(i).dstation) = RoadNetwork.Starters(1,car(i).dstation) + 1;
                        end
                    elseif car(i).state == IDLE
                        RoadNetwork.Starters(1,car(i).ostation) = RoadNetwork.Starters(1,car(i).ostation) + 1;
                    elseif car(i).state == DRIVING_TO_PICKUP
                        distToDest = norm(car(i).dpos - car(i).pos,2);
                        eta1 = distToDest / car(i).speedfactor;
                        cus = customer(car(i).passId);
                        if cus.ostation == cus.dstation
                            eta2 = 1;
                        else
                            eta2 = LinkTime(cus.ostation, cus.dstation);
                        end
                        tin = ceil(eta / (predictionStep * dt));
                        if tin < horizon && tin > 0
                            RoadNetwork.Starters(tin,car(i).dstation) = RoadNetwork.Starters(tin,car(i).dstation) + 1;
                        elseif tin < horizon
                            RoadNetwork.Starters(1,car(i).dstation) = RoadNetwork.Starters(1,car(i).dstation) + 1;
                        end
                    end
                end
                % load waiting customers into FlowsOut
                for st=1:numStations
                    for cu=1:length(station(st).custUnassigned)
                        custInd = station(st).custUnassigned(cu);
                        tout = 1;
                        Passengers.FlowsOut{tout}(customer(custInd).onode,customer(custInd).dnode) = Passengers.FlowsOut{tout}(customer(custInd).onode,customer(custInd).dnode) + 1; 
                    end
                end

                % run!
                [rebalanceQueue] = MPC_MCF(RoadNetwork,RebWeight,Passengers,Flags);
                for i=1:length(rebalanceQueue)
                    numRebTasks(t,i) = numRebTasks(t,i) + length(rebalanceQueue{i});
                end
            end
        else
            % calculate difference between the floor value and the actual
            % number
            vdesiredDifference = (sum(vown) - totalCustomers) - numStations*vdesired(1);
            [~, sortedIndex] = sort(vown);
            for i = 1:vdesiredDifference
                vdesired(sortedIndex(i)) = vdesired(sortedIndex(i)) + 1;
            end
            % form the sets S and T (sources and sinks) and the flows in/out
            sources = [];
            sinks = [];
            flowOut = [];
            flowIn = [];
            for i = 1:numStations
                tmpFlow = vexcess(i) - vdesired(i);
                if tmpFlow > 0
                    for j = 1:tmpFlow
                        % source
                        sources = [sources, StationNodeID(i)];
                        flowOut = [flowOut, 1];
                    end
                elseif tmpFlow < 0
                    for j = 1:-tmpFlow
                        % sink
                        sinks = [sinks, StationNodeID(i)];
                        flowIn = [flowIn, 1];
                    end
                end
            end
            %sum(flowOut)
            % find the link capacity left
            LinkCapacityLeft = max(RoadCap - LinkNumVehicles, 0);

            % call the single-commodity flow solver
            RebOutput = TIMulticommodityFlow(RoadGraph, LinkCapacityLeft, LinkTime , sources, sinks, flowIn, flowOut, 0, 0, 1);
            
            if PLOTREBFLAG
                
                TIPlotRebFlows(1,N,RoadGraph,NodesLocation,RebOutput,sources,sinks,0,1)
                
                
            end

            %This modifies the sources and sinks so the path planner knows that
            %some paths were modified
            S=length(sources);
            M=1;
            SoRelaxkl=@(k,l) N*N*M + N*N + S*(k-1) + l;
            SiRelaxkl=@(k,l) N*N*M + N*N + M*S + S*(k-1) + l;
            relFlowIn = flowIn -RebOutput(SoRelaxkl(1,1):SoRelaxkl(1,S))';
            relFlowOut= flowOut-RebOutput(SiRelaxkl(1,1):SiRelaxkl(1,S))';
            numVehiclesNotRebalancing(t) = sum(RebOutput(SoRelaxkl(1,1):SoRelaxkl(1,S)));
            rebPaths = TIRebPathDecomposition(RebOutput, length(RoadGraph), 1, sources, sinks, relFlowIn, relFlowOut);

            disp('how much didnt rebalance')
            sum(flowOut) - length(rebPaths)
        end
    else
        rebCount = rebCount + 1;        
    end
    
    % assign rebalancing vehicles
    if NAIVEFLAG
        for i = 1:numStations
            % rebalance cars as long as the station has cars available
            while ~isempty(rebalanceQueue{i}) && (~isempty(station(i).carIdle) || ~isempty(station(i).carOnRoad))
                currDest = rebalanceQueue{i}(1);
                % remove current destination from rebalanceQueue
                rebalanceQueue{i} = rebalanceQueue{i}(2:end);
                distToStation = norm(StationLocation(i,:) - StationLocation(currDest,:));
                if ~isempty(station(i).carIdle)
                    assignedCarID = station(i).carIdle(1);
                    shortestDistance = distToStation;
                else
                    shortestDistance = inf;
                end
                carOnRoadIndex = 0;
                if ~isempty(station(i).carOnRoad)
                    for j = 1:length(station(i).carOnRoad)
                        distToCar = norm(StationLocation(currDest,:) - car(station(i).carOnRoad(j)).pos);
                        if distToCar < shortestDistance
                            shortestDistance = distToCar;
                            assignedCarID = station(i).carOnRoad(j);
                            carOnRoadIndex = j;
                        end
                    end
                end
                if shortestDistance == inf
                    'ERROR: inf shortest distance in rebalancing'
                end
                % remove assigned car from station
                if ~isempty(station(i).carIdle) && assignedCarID == station(i).carIdle(1)
                    station(i).carIdle = station(i).carIdle(2:end);
                else
                    station(i).carOnRoad = [station(i).carOnRoad(1:carOnRoadIndex-1), station(i).carOnRoad(carOnRoadIndex+1:end)];
                end
                % remove car from station
                car(assignedCarID).dstation = currDest;
                car(assignedCarID).ostation = i;
                car(assignedCarID).state = REBALANCING;
                car(assignedCarID).dpos = StationLocation(currDest,:);
                tmpPath = car(assignedCarID).path;
                car(assignedCarID).path = findRoute(tmpPath(1), StationNodeID(currDest), LinkTime); 
                
                if length(tmpPath) > 1 && length(car(assignedCarID).path) > 1
                    if car(assignedCarID).path(2) ~= tmpPath(2)
                        % go a different way - make U-turn
                        car(assignedCarID).path = [tmpPath(2) car(assignedCarID).path];
                    end
                end
                if length(car(assignedCarID).path) > 1
                    LinkNumVehicles(car(assignedCarID).path(1), car(assignedCarID).path(2)) = LinkNumVehicles(car(assignedCarID).path(1), car(assignedCarID).path(2))+1;
                end
                % direction vector
                if length(car(assignedCarID).path) > 1 
                    tmpNode2 = car(assignedCarID).path(2);
                    tmpNode1 = car(assignedCarID).path(1);
                    car(assignedCarID).direction = (NodesLocation(tmpNode2,:) - NodesLocation(tmpNode1,:)); car(assignedCarID).direction=car(assignedCarID).direction/norm(car(assignedCarID).direction);
                end
                % destination station add this car
                station(currDest).carOnRoad = [station(currDest).carOnRoad, assignedCarID];
            end
        end

    else
        % give the rebPaths to the vehicles
        nextRebPaths = cell(0);
        c_reb = 1;
        for i = 1:length(rebPaths)
            if isempty(rebPaths{i})
                continue
            end
            orig = rebPaths{i}{1,1}(1,1);
            stationId = find(StationNodeID == orig,1);
            %dest = rebPaths{1}(end);
            if ~isempty(station(stationId).carIdle)
                % assign car to rebalance path
                tmpCar = station(stationId).carIdle(1);
                car(tmpCar).state = REBALANCING;
                car(tmpCar).path = rebPaths{i}{1,1}(:,1)';
                car(tmpCar).dpos = NodesLocation(rebPaths{i}{1,1}(end,1),:);
                car(tmpCar).direction = (NodesLocation(car(tmpCar).path(2),:) - NodesLocation(car(tmpCar).path(1),:)); car(tmpCar).direction=car(tmpCar).direction/norm(car(tmpCar).direction);
                car(tmpCar).dstation = find(StationNodeID == car(tmpCar).path(end),1);
                station(stationId).carIdle = station(stationId).carIdle(2:end);
                LinkNumVehicles(car(tmpCar).path(1), car(tmpCar).path(2)) = LinkNumVehicles(car(tmpCar).path(1), car(tmpCar).path(2))+1;
                station(car(tmpCar).dstation).carOnRoad = [station(car(tmpCar).dstation).carOnRoad tmpCar];
            else
                % save this for next time
                nextRebPaths{c_reb} = rebPaths{i};
                c_reb = c_reb + 1;
            end

        end

        rebPaths = nextRebPaths;
    
    end
    
    % move vehicles
    for i = 1:v
        if car(i).state == SELF_LOOP
            car(i).time_left=car(i).time_left-dt;
        else
            if length(car(i).path) > 1
                if norm(NodesLocation(car(i).path(2),:)-car(i).pos) < car(i).speedfactor * LinkSpeed(car(i).path(1), car(i).path(2))*dt
                    % need to go to the next node (or stop)
                    
                    residualDistance = car(i).speedfactor * LinkSpeed(car(i).path(1), car(i).path(2))*dt - norm(NodesLocation(car(i).path(2),:)-car(i).pos);
                    
                    if length(car(i).path) > 2
                        LinkNumVehicles(car(i).path(1), car(i).path(2)) = LinkNumVehicles(car(i).path(1), car(i).path(2)) - 1;
                        car(i).path = car(i).path(2:end);
                        car(i).direction = (NodesLocation(car(i).path(2),:) - NodesLocation(car(i).path(1),:)); car(i).direction=car(i).direction/norm(car(i).direction);
                        car(i).pos = NodesLocation(car(i).path(1),:) + car(i).direction*residualDistance;
                        if isnan(car(i).pos)
                            'ERROR: car position is NaN'
                        end
                        LinkNumVehicles(car(i).path(1), car(i).path(2)) = LinkNumVehicles(car(i).path(1), car(i).path(2)) + 1;
                    else
                        car(i).pos = NodesLocation(car(i).path(2),:);
                        if isnan(car(i).pos)
                            'Error: car position is NaN'
                        end
                    end
                    
                else
                    % don't need to go to the next node
                    car(i).pos = car(i).pos + car(i).direction * car(i).speedfactor * LinkSpeed(car(i).path(1), car(i).path(2))*dt;
                    if isnan(car(i).pos)
                        'Error: car position is NaN'
                    end
                end
                
            end
        end
    end
    
    
    % update customer waiting/service times
    for i = 1:cc-1
        if customer(i).pickedup == 0
            customer(i).waitTime = customer(i).waitTime + dt;
        elseif customer(i).pickedup == 1 && customer(i).delivered == 0
            customer(i).serviceTime = customer(i).serviceTime + dt;
        end
    end
    
    % collect customer travel time data
    for i = 1:numStations
        for j = 1:length(station(i).custId)
            if customer(station(i).custId(j)).pickedup == 0
                avgWaitTimes(t) = avgWaitTimes(t) + customer(station(i).custId(j)).waitTime;
                cumNumCustomers(t) = cumNumCustomers(t) + 1;
            end
        end
    end
    if cumNumCustomers(t) > 0
        avgWaitTimes(t) = avgWaitTimes(t)./cumNumCustomers(t);
    end
    
    % collect number of vehicles doing stuff
    for i = 1:v
        if car(i).state ~= IDLE
            numVehiclesBusy(t) = numVehiclesBusy(t) + 1;
            if car(i).state == REBALANCING
                numVehiclesRebalancing(t) = numVehiclesRebalancing(t) + 1;
            end
            if car(i).state == DRIVING_TO_PICKUP
                numVehiclesDrivingToPickup(t) = numVehiclesDrivingToPickup(t)+1;
            end
            if car(i).state == DRIVING_TO_DEST
                numVehiclesDrivingToDest(t) = numVehiclesDrivingToDest(t)+1;
            end
            if car(i).state == DRIVING_TO_STATION
                numVehiclesDrivingToStation(t) = numVehiclesDrivingToStation(t)+1;
            end
            if car(i).state == SELF_LOOP
                numVehiclesOnSelfLoop(t) = numVehiclesOnSelfLoop(t)+1;
            end
        else
            numVehiclesIdle(t) = numVehiclesIdle(t)+1;
        end
    end
end

if PLOTFLAG
    close(writerObj);
end

%%
numDelivered = 0;
totalServiceTime = 0;
totalWaitTime = 0;
for i = 1:cc-1
    if customer(i).delivered
        numDelivered = numDelivered + 1;
        totalServiceTime = totalServiceTime + customer(i).serviceTime;
        totalWaitTime = totalWaitTime + customer(i).waitTime;
    end
end
meanServiceTime = totalServiceTime/numDelivered;

% histogram of waiting times
allWaitTimes = zeros(numDelivered, 1);
ccc = 1;
for i = 1:cc-1
    if customer(i).delivered
        allWaitTimes(ccc) = customer(i).waitTime;
        ccc = ccc+1;
    end
end

allServiceTimes = zeros(numDelivered, 1);
ccc = 1;
for i = 1:cc-1
    if customer(i).delivered
        allServiceTimes(ccc) = customer(i).serviceTime;
        ccc = ccc+1;
    end
end


savename = [SIMNAME,num2str(v),'v_rev3'];
save(savename);
function FlowsOut = predictDemand(currentTime, horizon, numStations, predictor)
	FlowsOut = cell(1, horizon);
	for t=1:horizon
		FlowsOut{t} = sparse(numStations,numStations);
		if t <= length(predictor{currentTime})
			% there is a prediction for that time
			% predictor here is a cell of cells
			FlowsOut{t} = FlowsOut{t} + predictor{currentTime}{t};
		end
	end
end
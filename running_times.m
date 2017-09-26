horizons = [3, 6, 12, 24, 48];
% look_back, look_forward, mean_wait, median_wait
times = [zeros(length(horizons)*length(horizons), 4);]
r = 1;
for i=1:length(horizons)
	for j=1:length(horizons)
		results(r,:) = [horizons(i) horizons(j) 0 0 ];
		r = r + 1;
	end
end

times = [];

for r=1:25
	ind = length(horizons)*length(horizons) + 1 - r
	hors = struct;
	a = results(ind,:);
	suffix = strcat(num2str(a(1)),'b_',num2str(a(2)),'f_np')
	FNAME = strcat('ignored_assets/','MPC-LSTM-MILP-', suffix,'5000v_rev3.mat');
	if exist(FNAME,'file') == 2
		milpOutputs = load(FNAME,'milpOutputs');
		milpOutputs = milpOutputs.milpOutputs;
		mask = milpOutputs(:,1) > 0;
		times = [times ; milpOutputs(mask,1)];
		fprintf('FOUND FILE!!')
	end
end

fprintf('Number of samples: %d \n',length(times))
fprintf('Mean: %d \n',mean(times))
fprintf('Median: %d \n',median(times))
fprintf('STD: %d \n',std(times))
fprintf('Max: %d \n',max(times))

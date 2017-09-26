horizons = [3, 6, 12, 24, 48];
% look_back, look_forward, mean_wait, median_wait
results = zeros(length(horizons)*length(horizons), 4);
r = 1;
for i=1:length(horizons)
	for j=1:length(horizons)
		results(r,:) = [horizons(i) horizons(j) 0 0 ];
		r = r + 1;
	end
end

minind = 1;
maxind = 25;
for r=minind:maxind
	ind = r%maxind + 1 - r
	hors = struct;
	a = results(ind,:);
	suffix = strcat(num2str(a(1)),'b_',num2str(a(2)),'f_np')
	FNAME = strcat('ignored_assets/','MPC-LSTM-MILP-', suffix,'5000v_rev3.mat');
	if exist(FNAME,'file') == 2
		res = load(FNAME,'mean_wait', 'median_wait');
		mean_wait = res.mean_wait;
		median_wait = res.median_wait;
		fprintf('FOUND FILE!!')
	else
		hors.lb = a(1);
		hors.lf = a(2);
		[mean_wait, median_wait] = run_simulation(hors);
	end
	results(ind,:) = [a(1) a(2) mean_wait median_wait];
end

%r = 1;
%for i=1:length(horizons)
%	parfor j=1:length(horizons)
%		hors = struct;
%		hors.lb = num2str(horizons(i));
%		hors.lf = num2str(horizons(j));
%		[mean_wait, median_wait] = run_simulation(hors);
%		results(r,:) = [horizons(i), horizons(j), mean_wait, median_wait];
%		r = r + 1;
%	end
%end

save('explore_horizons_results.mat', 'results');
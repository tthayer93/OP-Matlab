function [time_data] = SCOP_to_CMDP_adaptive_time_data(vertex_clusters, t_max, num_trials, adaptive_type)
%SCOP_TO_CMDP_ADAPTIVE_TIME_DATA Simulates time data for following a specified route with stochastic costs
%
%	Version: 1.0
%	Date: 04/28/20
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function produces simulated arrival time data for following an input route with stochastic costs, which is used to solve the Stochastic Cost Orienteering Problem (SCOP) by computing a Constrained Markov Decision Proccess (CMDP) policy where the time discretization is chosen with respect to the simulated time data. First used in "A Resolution Adaptive Algorithm for the Stochastic Orienteering Problem with Chance Constraints" by Thayer, Carpin
%	Assumptions:
%		The clustered SCOP starts at clustered vertex 1 at time 0.
%		The goal is at clustered vertex N.
%		The failure states are encapsulated by vertices N+1 to N+N, where vertex N+i represents failing when leaving clustered vertex i.
%		The absorbing state is encapsulated by vertex N+N+2.
%		The failure state is achieved when time > t_max.
%		Clustered vertex i can only travel to clustered vertex j such that i < j <= N.
%		The costs for traveling from clustered vertex x to clustered vertex y is defined by the associated minimum cost and combined distribution.
%	Inputs:
%		vertex_clusters: an upper triangular cell matrix of size NxNx3 that contains the following information
%			cell (x,y,1): a vector of the sequence of verticies to traverse from the end of the current cluster x (row) to the end of cluster y (column)
%			cell (x,y,2): a sequence (whos lengths match the corresponding sequence of verticies) of cells each containing
%				1: the minimum cost for the associated vertex transition
%				2: the type of single parameter distribution used to model the cost for the associated vertex transition
%				3: the parameter of the distribution for the associated vertex transition
%			cell (x,y,3): a cell containing the combined distribution of the sequence of verticies for the associated vertex cluster transition with
%				1: the minimum cost for the associated vertex cluster transition
%				2: the type of distribution used to model the cost for the associated vertex cluster transition
%				3+: the parameters of the distribution for the associated vertex cluster transition
%		t_max: the time budget of the SCOP
%		num_trials: the number of time data points to simulate for each vertex_cluster
%		adaptive_type:
%	Outputs:
%		time_data: a cell vector of length N where each cell contains the simulated arrival time data for the corresponding vertex_cluster, including an additional num_trials/2 equally spaced times from 0 to t_max

	%% Check inputs
	if exist('adaptive_type') && ~isempty(adaptive_type)
		if adaptive_type < 0
			%warning('adaptive_type must be positive.');
			time_data = [];
			return;
		end
	else
		adaptive_type = 0;
	end
	%% Determine policy with adaptive time
	time_data = cell(size(vertex_clusters, 1), 1);
	time_data{1} = zeros(num_trials, 1);
	for i=2:size(vertex_clusters, 1)
		j = 1;
		this_dist = vertex_clusters{i-j, i, 3};
		while isempty(this_dist)
		%while this_dist{3} == 0
			j = j + 1;
			this_dist = vertex_clusters{i-j, i, 3};
		end
		idx = randperm(size(time_data{i-1}, 1), num_trials);
		time_data{i} = time_data{i-1}(idx) + this_dist{1} + random(this_dist{2:end}, [num_trials, 1]);
	end
	for i=1:size(vertex_clusters, 1)
		if adaptive_type > 0
			time_data{i} = [time_data{i}; [0:(t_max/(num_trials/adaptive_type)):t_max]'];
		end
		time_data{i} = time_data{i}((time_data{i} <= t_max) & (time_data{i} >= 0));
	end
	
end
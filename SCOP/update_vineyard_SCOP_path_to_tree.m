function [vertex_cluster_rewards, vertex_clusters] = update_vineyard_SCOP_path_to_tree(vertex_cluster_rewards, vertex_clusters, t_max, tour, states, state_action, rho, reward_map, row_distance, vine_distance, start_vertex, end_vertex, cost_dist_type, cost_min_frac, heuristic, this_size_clusters)
%UPDATE_VINEYARD_SCOP_PATH_TO_TREE Creates a SCOP path tree from a previously solved SCOP on a vineyard graph
%
%	Version: 1.0
%	Date: 12/02/20
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function takes a vineyard Stochastic Cost Orienteering Problem (SCOP) path that has been solved using a Constrained Markov Decision Process (CMDP) and converts it to a SCOP path tree, using a heuristic to reduce the number of branches added. A vineyard adaptation of the algorithm in "An Adaptive Method for the Stochastic Orienteering Problem" by Thayer, Carpin.
%	Assumptions:
%		The clustered SCOP starts at clustered vertex 1 at time 0.
%		The goal is at clustered vertex N.
%		Clustered vertex i can only travel to clustered vertex j such that i < j <= N and i is on the same branch as or the originating branch of j.
%		The costs for traveling from clustered vertex x to clustered vertex y is defined by the associated minimum cost and combined distribution.
%	Inputs:
%		vertex_cluster_rewards: a vector containing the reward (positive values) for visiting each clustered vertex
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
%		tour: the initial OP tour on which the SCOP is based
%		states: a table containing the relevant information for each state, with the following format
%			state | clustered_vertex | start_time | end_time
%		state_action: A matrix containing all possible state-maneuver pairs, with the following format
%			state | maneuver
%		rho: a vector containing the computed rho values of the CMDP
%		reward_map: A matrix of size num_rows*num_vines_per_row containing the reward of each vertex
%		row_distance: The distance between each row of vines
%		vine_distance: The distance between each vine in the rows
%		start_vertex: the first vertex in the route
%		end_vertex: the last vertex in the route
%		cost_dist_type: the type of single parameter distribution used to model costs between adjacent vertices in the vineyard
%		cost_min_frac: the fraction of cost for the minimum bound of each possible transition, such that cost*cost_min_frac + distribution(1-cost_min_frac) = cost of transition
%		heuristic: An integer value from 0 to inf that determines how the branching heuristic should work
%			heuristic = 0: No branches are added and the function immediately returns the original SCOP
%			heuristic = 1: All states where a policy action exists are used to determine branches
%			heuristic = 2: All state_action pairs with a rho value larger than the mean are used to determine branches
%			heuristic >= 3: The input value is utilized such that the top k values of rho are used to determine branches
%		this_size_clusters: the maximum number of vertex clusters to group together into a single vertex cluster
%	Outputs:
%		vertex_cluster_rewards: a vector containing the reward (positive values) for visiting each clustered vertex
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

	old_vertex_clusters = vertex_clusters;
	old_vertex_cluster_rewards = vertex_cluster_rewards;
	if heuristic == 0
		return
	elseif heuristic == 1
		policy_null = compute_CMDP_policy(rho, state_action);
		temp = state_action(policy_null(:, 1) > 0, :);
		temp = temp(temp(:, 2) <= length(tour), :);
		states_table = table2array(states);
		jump_states = temp(logical(temp(:, 2) - states_table(temp(:, 1), 2) - 1), 1);
		jump_states = jump_states(jump_states > 1);
	elseif heuristic == 2
		min_rho = mean(rho);
		temp = state_action(rho > min_rho, :);
		temp = temp(temp(:, 2) <= length(tour), :);
		states_table = table2array(states);
		jump_states = temp(logical(temp(:, 2) - states_table(temp(:, 1), 2) - 1), 1);
		jump_states = jump_states(jump_states > 1);
	elseif heuristic >= 3
		policy_null = compute_CMDP_policy(rho, state_action);
		sa_idx = (state_action(:, 2) <= length(tour));
		policy_null = policy_null(sa_idx);
		rho_null = rho(sa_idx);
		idx = (policy_null(:, 1) > 0);
		temp_rho = rho_null(idx);
		sa_temp = state_action(sa_idx, :);
		temp = sa_temp(idx, :);
		states_table = table2array(states);
		jump_idxs = logical(temp(:, 2) - states_table(temp(:, 1), 2) - 1) & (temp(:, 1) > 1);
		jump_states = temp(jump_idxs, 1);
		jump_rho = temp_rho(jump_idxs);
		jump_rho(states_table(jump_states, 2) >= size(old_vertex_clusters, 1)) = 0;
		if sum(logical(jump_rho)) < heuristic
			jump_states = jump_states(logical(jump_rho), :);
			jump_rho = jump_rho(logical(jump_rho), :);
		end
		[kept ,temp_idxs] = max_k(jump_rho, heuristic);
		jump_states = jump_states(temp_idxs, :);
	end
	for i=1:length(jump_states)
		this_start_cluster = states_table(jump_states(i), 2);
		this_start_vertex = old_vertex_clusters{1, this_start_cluster, 1}(end);
		this_budget = t_max - states_table(jump_states(i), 4);
		removable_vertices = tour(1:(find(tour == this_start_vertex) - 1));
		this_rewards = reward_map;
		this_rewards(floor((removable_vertices - 1) / size(reward_map, 2)) + 1, mod(removable_vertices - 1, size(reward_map, 2)) + 1) = 0;
		[this_total_cost, this_total_reward, this_tour] = greedy_partial_row_avoidance(vine_distance, row_distance, this_rewards, this_start_vertex, end_vertex, this_budget);
		if isempty(this_tour)
			continue;
		end
		temp = [start_vertex];
		for j=1:this_start_cluster-1
			temp = [temp; old_vertex_clusters{j, j+1}];
		end
		new_tour = [temp; this_tour(2:end)];
		if (length(new_tour) == length(tour)) && all(new_tour == tour)
			continue;
		end
		split_begin = 1;
		while new_tour(split_begin) == tour(split_begin)
			split_begin = split_begin + 1;
		end
		split_end = 0;
		while new_tour(end-split_end) == tour(end-split_end)
			split_end = split_end + 1;
		end
		split_end = split_end + 1;
		if length(new_tour)-split_end >= split_begin
			%cluster vineyard rows and such
			[this_vertex_cluster_rewards, this_vertex_clusters] = vineyard_route_to_SCOP(vine_distance, row_distance, this_rewards, new_tour, cost_dist_type, cost_min_frac);
			[this_vertex_cluster_rewards, this_vertex_clusters] = cluster_SCOP(this_vertex_cluster_rewards, this_vertex_clusters, this_size_clusters);
			temp_length = 0;
			for j=1:size(this_vertex_clusters)
				temp_length = temp_length + length(this_vertex_clusters{j, j+1, 1});
				if temp_length > split_begin
					split_begin = j;
					break;
				end
			end
			new_vertex_cluster_rewards = [vertex_cluster_rewards(1:end-1), this_vertex_cluster_rewards(split_begin:end)];
			new_vertex_clusters = cell(length(new_vertex_cluster_rewards), length(new_vertex_cluster_rewards), 3);
			new_vertex_clusters(1:length(vertex_cluster_rewards)-1, 1:length(vertex_cluster_rewards)-1, :) = vertex_clusters(1:end-1, 1:end-1, :);
			new_vertex_clusters(1:length(vertex_cluster_rewards)-1, end, :) = vertex_clusters(1:end-1, end, :);
			new_vertex_clusters(1:(split_begin-1), length(vertex_cluster_rewards):end-1, :) = this_vertex_clusters(1:(split_begin-1), split_begin:end-1, :);
			new_vertex_clusters(length(vertex_cluster_rewards):end, length(vertex_cluster_rewards):end, :) = this_vertex_clusters(split_begin:end, split_begin:end, :);
			vertex_clusters = new_vertex_clusters;
			vertex_cluster_rewards = new_vertex_cluster_rewards;
		end
	end
% 	for i=1:length(vertex_cluster_rewards)
% 		for j=(i+1):length(vertex_cluster_rewards)
% 			if size(vertex_clusters{i, j, 3}, 1) == 0
% 				vertex_clusters{i, j, 2} = [num2cell(0), {'None'}, 0];
% 				vertex_clusters{i, j, 3} = {0, 'Gamma', 0, 0};
% 			end
% 		end
% 	end
	
end


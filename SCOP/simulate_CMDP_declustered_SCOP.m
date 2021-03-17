function [tour, state_path, maneuvers, total_reward, failure] = simulate_CMDP_declustered_SCOP(policy, state_action, state_transition_table, states, vertex_clusters, vertex_rewards, start_vertex, failure_penalties)
%SIMULATE_CMDP_DECLUSTERED_SCOP Simulates a run of a precomputed CMDP policy for a clustered SCOP in the original declustered space
%
%	Version: 1.0
%	Date: 04/11/19
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function simulates a single run of a precomputed Constrained Markov Decision Process (CMDP) policy for a particular Stochastic Cost Orienteering Problem (SCOP) that has been clustered in the original declustered space.
%	Assumptions:
%		The SCOP starts at vertex cluster 1 at time 0.
%		The goal is at vertex cluster N.
%		The failure state is encapsulated by vertex cluster N+1.
%		The absorbing state is encapsulated by vertex cluster N+2.
%		The failure state is achieved when time > t_max.
%		Vertex cluster i can only travel to vertex cluster j such that i < j <= N+1.
%		The costs for traveling from vertex cluster x to vertex cluster y is defined by the cost to travel from the last vertex in cluster x to the first vertex in cluster y and follow the sequence of verticies in cluster y.
%	Inputs:
%		policy: a vector containing the optimal policy, which follows the corrresponding state_action
%		state_action: A matrix containing all possible state-maneuver pairs, with the following format
%			state | maneuver
%		state_transition_table: a table containing all possible state transitions (one per row) and associated values where V >= 0, with the following format
%			state | maneuver | next_state | probability | value_0 | ... | value_V
%		states: a table containing the relevant information for each state, with the following format
%			state | vertex | start_time | end_time
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
%		vertex_rewards: a table containing the reward for visiting each vertex, with the following format
%			vertex | reward
%		start_vertex: the vertex associated with vertex cluster 1, at which the tour should start
%		failure_penalties: a vector containing the reward (negative values) for failing when leaving each clustered vertex
%	Outputs:
%		tour: a vector containing the sequence of vertices visited by the simulation
%		state_path: a vector containing the sequence of states visited by the simulation
%		maneuvers: a vector containing the sequence of maneuvers chosen by the policy in simulation
%		total_reward: the reward accumulated by the simulation
%		failure: a binary variable that is true if the simulation reaches the failure state

	%% Check inputs
	if ~exist('failure_penalties') || isempty(failure_penalties)
		failure_penalties = zeros(size(vertex_clusters, 1), 1);
	else
		if length(failure_penalties) == 1
			failure_penalties = failure_penalties * ones(size(vertex_clusters, 1), 1);
		elseif length(failure_penalties) ~= size(vertex_clusters, 1)
			error('Length of failure_penalties must equal 0, 1, or N.');
		end
	end
	%% Simulate the policy
	tour = [start_vertex];
	state_transition_table = table2array(state_transition_table);
	states = table2array(states);
	maneuvers = zeros(length(policy), 1);
	maneuvers_counter = 0;
	state_path = zeros(length(policy), 1);
	state_path(1) = states(1, 1);
	state_path_counter = 1;
	time = 0;
	this_clustered_vertex = 1;
	failure = 0;
	failure_penalty = 0;
	while ~any(state_path(state_path_counter) == states(end, 1))
		%% Expand data structures if needed
		if maneuvers_counter == length(maneuvers)
			maneuvers = [maneuvers; zeros(length(policy), 1)];
		end
		if state_path_counter == length(state_path)
			state_path = [state_path; zeros(length(policy), 1)];
		end
		%% Find next maneuver and possible transitions
		idx = (state_action(:, 1) == state_path(state_path_counter)) & (policy > 0);
		if sum(idx) > 1
			maneuvers(maneuvers_counter + 1) = randsample(state_action(idx, 2), 1, 1, policy(idx));
			maneuvers_counter = maneuvers_counter + 1;
		elseif sum(idx) == 1
			maneuvers(maneuvers_counter + 1) = state_action(idx, 2);
			maneuvers_counter = maneuvers_counter + 1;
		else
			error('No policy exists for the current state. Try increasing the number of time steps');
		end
		temp_table = state_transition_table((state_transition_table(:, 1) == state_path(state_path_counter)) & (state_transition_table(:, 2) == maneuvers(maneuvers_counter)), :);
		%% Get next state
		%next_clustered_vertex = states(states(:, 1) == min(temp_table(:, 3)), 2);
		next_clustered_vertexes = states(any(states(:, 1) == temp_table(:, 3)', 2), 2);
		next_clustered_vertex = next_clustered_vertexes(1);
		fail_vertex = next_clustered_vertexes(end);
		if next_clustered_vertex <= size(vertex_clusters, 1)
			[this_path, this_dist] = vertex_clusters{this_clustered_vertex, next_clustered_vertex, 1:2};
% 			for i=1:length(this_path)
% 				time = time + this_dist{i, 1} + random(this_dist{i, 2:end});
% 				if time > states(end-1, 3)
% 					i = i - 1;
% 					break;
% 				end
% 			end
% 			tour = [tour; this_path(1:i)];
			%
			time_vec = cumsum(vertcat(this_dist{:, 1}) + random('exponential', vertcat(this_dist{:, 3})));
			i = find(time + time_vec > states(end-1, 3), 1);
			if isempty(i)
				tour = [tour; this_path];
				time = time + time_vec(end);
			else
				tour = [tour; this_path(1:(i-1))];
				time = time + time_vec(i);
			end
			%
			idx = ((states(:, 2) == next_clustered_vertex) & (states(:, 3) < time) & (states(:, 4) >= time));
			if sum(idx) == 0
				%idx = ((states(:, 3) < time) & (states(:, 4) >= time));
				idx = (states(:, 2) == fail_vertex);
				failure = 1;
				failure_penalty = failure_penalties(fail_vertex - size(vertex_clusters, 1));
			end
		else
			idx = (states(:, 2) == next_clustered_vertex);
		end
		state_path(state_path_counter + 1) = states(idx, 1);
		state_path_counter = state_path_counter + 1;
		this_clustered_vertex = states(idx, 2);
	end
	%% Shrink data structures to fit
	maneuvers((maneuvers_counter + 1):end) = [];
	state_path((state_path_counter + 1):end) = [];
	%% Add unclustered reward
	total_reward = sum(vertex_rewards(unique(tour), 2)) - failure_penalty;
	
end


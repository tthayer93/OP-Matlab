function [state_path, maneuvers, total_reward, failure, time_path] = simulate_CMDP_clustered_SCOP(policy, state_action, state_transition_table, states, vertex_cluster_rewards, vertex_clusters, failure_penalties)
%SIMULATE_CMDP_CLUSTERED_SCOP Simulates a run of a precomputed CMDP policy for a clustered SCOP
%
%	Version: 1.1
%	Date: 24/1/20
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function simulates a single run of a precomputed Constrained Markov Decision Process (CMDP) policy for a particular Stochastic Cost Orienteering Problem (SCOP).
%	Assumptions:
%		The SCOP starts at clustered vertex 1 at time 0.
%		The goal is at clustered vertex N.
%		The failure state is encapsulated by vertex N+1.
%		The absorbing state is encapsulated by vertex N+2.
%		The failure state is achieved when time > t_max.
%		Clustered vertex i can only travel to clustered vertex j such that i < j <= N+1.
%		The costs for traveling from vertex cluster x to vertex cluster y is defined by the cost to travel from the last vertex in cluster x to the first vertex in cluster y and follow the sequence of verticies in cluster y.
%	Inputs:
%		policy: a vector containing the optimal policy, which follows the corrresponding state_action
%		state_action: A matrix containing all possible state-maneuver pairs, with the following format
%			state | maneuver
%		state_transition_table: a table containing all possible state transitions (one per row) and associated values where V >= 0, with the following format
%			state | maneuver | next_state | probability | value_0 | ... | value_V
%		states: a table containing the relevant information for each state, with the following format
%			state | clustered vertex | start_time | end_time
%		vertex_cluster_rewards: a vector containing the reward for visiting each vertex cluster
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
%		failure_penalties: a vector containing the reward (negative values) for failing when leaving each clustered vertex
%	Outputs:
%		state_path: a vector containing the sequence of states visited by the simulation
%		maneuvers: a vector containing the sequence of maneuvers chosen by the policy in simulation
%		total_reward: the reward accumulated by the simulation
%		failure: a binary variable that is true if the simulation reaches the failure state
%		time_path: a vector containing the simulated time for reaching each state

	%% Check inputs
	if ~exist('failure_penalties') || isempty(failure_penalties)
		failure_penalties = zeros(length(vertex_cluster_rewards), 1);
	else
		if length(failure_penalties) == 1
			failure_penalties = failure_penalties * ones(length(vertex_cluster_rewards), 1);
		elseif length(failure_penalties) ~= length(vertex_cluster_rewards)
			error('Length of failure_penalties must equal 0, 1, or N.');
		end
	end
	%% Simulate the policy
	state_transition_table = table2array(state_transition_table);
	states = table2array(states);
	maneuvers = zeros(length(policy), 1);
	maneuvers_counter = 0;
	state_path = zeros(length(policy), 1);
	state_path(1) = states(1, 1);
	state_path_counter = 1;
	this_vertex = 1;
	total_reward = vertex_cluster_rewards(this_vertex);
	failure = 0;
	time_path = zeros(length(policy), 1);
	time_path(1) = 0;
	time_path_counter = 1;
	while ~any(state_path(state_path_counter) == states(end, 1))
		%% Expand data structures if needed
		if maneuvers_counter == length(maneuvers)
			maneuvers = [maneuvers; zeros(length(policy), 1)];
		end
		if state_path_counter == length(state_path)
			state_path = [state_path; zeros(length(policy), 1)];
		end
		if time_path_counter == length(time_path)
			time_path = [time_path; zeros(length(policy), 1)];
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
			error('No policy exists for the current state. Check transition probabilities or try increasing the number of time steps.');
		end
		temp_table = state_transition_table((state_transition_table(:, 1) == state_path(state_path_counter)) & (state_transition_table(:, 2) == maneuvers(maneuvers_counter)), :);
		%% Get next state
		%next_vertex = states(states(:, 1) == min(temp_table(:, 3)), 2);
		next_vertexes = states(any(states(:, 1) == temp_table(:, 3)', 2), 2);
		next_vertex = next_vertexes(1);
		fail_vertex = next_vertexes(end);
		if next_vertex <= size(vertex_clusters, 1)
			this_dist = vertex_clusters{this_vertex, next_vertex, 3};
			time_path(time_path_counter + 1) = time_path(time_path_counter) + this_dist{1} + random(this_dist{2:end});
			time_path_counter = time_path_counter + 1;
			idx = ((states(:, 2) == next_vertex) & (states(:, 3) < time_path(time_path_counter)) & (states(:, 4) >= time_path(time_path_counter)));
			if sum(idx) == 0
				%idx = ((states(:, 3) < time_path(time_path_counter)) & (states(:, 4) >= time_path(time_path_counter)));
				idx = (states(:, 2) == fail_vertex);
				failure = 1;
				total_reward = total_reward - failure_penalties(fail_vertex - size(vertex_clusters, 1));
			end
		else
			idx = (states(:, 2) == next_vertex);
			time_path(time_path_counter + 1) = time_path(time_path_counter);
			time_path_counter = time_path_counter + 1;
		end
		state_path(state_path_counter + 1) = states(idx, 1);
		state_path_counter = state_path_counter + 1;
		this_vertex = states(idx, 2);
		%% Add reward
		if this_vertex <= length(vertex_cluster_rewards)
			total_reward = total_reward + vertex_cluster_rewards(this_vertex);
		end
	end
	%% Shrink data structures to fit
	maneuvers((maneuvers_counter + 1):end) = [];
	state_path((state_path_counter + 1):end) = [];
	time_path((time_path_counter + 1):end) = [];
	
end


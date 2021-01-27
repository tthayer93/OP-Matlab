function [states, state_transition_table, initial_state_distribution, absorbing_state] = clustered_SCOP_to_CMDP_multifail(vertex_cluster_rewards, vertex_clusters, t_max, num_time_steps, num_state_steps, failure_penalties)
%CLUSTERED_SCOP_TO_CMDP Creates a CMDP state_transition_table from a clustered SCOP
%
%	Version: 1.0
%	Date: 05/11/19
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function creates a Constrained Markov Decision Process (CMDP) formulation from a clustered Stochastic Cost Orienteering Problem (SCOP) that has been solved in expectation as an ordinary Orienteering Problem.
%	Assumptions:
%		The clustered SCOP starts at clustered vertex 1 at time 0.
%		The goal is at clustered vertex N.
%		The failure states are encapsulated by vertices N+1 to N+N, where vertex N+i represents failing when leaving clustered vertex i.
%		The absorbing state is encapsulated by vertex N+N+2.
%		The failure state is achieved when time > t_max.
%		Clustered vertex i can only travel to clustered vertex j such that i < j <= N.
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
%		num_time_steps: the number of time steps to discretize time with
%		num_state_steps: the number of steps used to discretize each state when computing transition probabilities to future states
%		failure_penalties: a vector containing the reward (negative values) for failing when leaving each clustered vertex
%	Outputs:
%		states: a table containing the relevant information for each state, with the following format
%			state | clustered_vertex | start_time | end_time
%		state_transition_table: a table containing all possible state transitions (one per row) and associated values (value_0 = -1 * reward for leaving a clustered vertex or -1 * penalty for leaving a failure vertex, value_1 = 1 for leaving the failure vertex and 0 otherwise), with the following format
%			state | maneuver | next_state | probability | value_0 | value_1
%		initial_state_distribution: a tuple containing the initial state and probability, i.e. (1,1)
%		absorbing_states: an integer denoting the absorbing state
	
	%% Check inputs
	if any(length(vertex_cluster_rewards) ~= size(vertex_clusters, 1))
		error('Number of vertices must be consistent for all inputs.');
	end
	if ~exist('num_state_steps')
		num_state_steps = 10;
	elseif isempty(num_state_steps)
		num_state_steps = 10;
	elseif num_state_steps <= 1
		num_state_steps = 2;
	end
	if ~exist('failure_penalties') || isempty(failure_penalties)
		failure_penalties = zeros(length(vertex_cluster_rewards), 1);
	else
		if length(failure_penalties) == 1
			failure_penalties = failure_penalties * ones(length(vertex_cluster_rewards), 1);
		elseif length(failure_penalties) ~= length(vertex_cluster_rewards)
			error('Length of failure_penalties must equal 0, 1, or N.');
		end
	end
	%% Create states
	num_vertex = length(vertex_cluster_rewards);
	failure_vertexes = [(num_vertex+1):(num_vertex+num_vertex)]';
	absorbing_vertex = num_vertex + num_vertex + 1;
	states = zeros((num_vertex-1)*num_time_steps + num_vertex + 2, 4);
	states(1, :) = [1, 1, 0, 0];
	counter = 1;
	time_step = t_max/num_time_steps;
	for i=2:num_vertex
		for j=0:time_step:t_max-time_step
			counter = counter + 1;
			states(counter, :) = [counter, i, j, j+time_step];
		end
	end
	failure_states = [(counter+1):(counter+num_vertex)]';
	states((counter+1):(counter+num_vertex), :) = [failure_states, failure_vertexes, t_max * ones(num_vertex, 1), inf * ones(num_vertex, 1)];
	absorbing_state = counter + num_vertex + 1;
	states(counter + num_vertex + 1, :) = [absorbing_state, absorbing_vertex, inf, inf];
	%% Create initial state distribution
	initial_state_distribution = [1, 1];
	%% Create state transition table
	state_transition_table = zeros((num_vertex-1)*num_time_steps, 6);
	counter = 0;
	for i = 1:size(states, 1)-2
		this_vertex = states(i, 2);
		if this_vertex < num_vertex
			this_time = states(i, 3);
			temp_table = states((states(:, 2) > this_vertex) & all(states(:, 2) < failure_vertexes', 2), :);
			%temp_table = temp_table(temp_table(:, 3) >= this_time, :);
			temp_table = temp_table(temp_table(:, 4) >= this_time, :);
			next_vertexes = unique(temp_table(:, 2));
			for j=1:length(next_vertexes)
				temp2_table = temp_table(temp_table(:, 2) == next_vertexes(j), :);
				%% Expand data structure if needed
				if (counter + size(temp2_table, 1) + 2) > size(state_transition_table, 1)
					state_transition_table = [state_transition_table; zeros((num_vertex-1)*num_time_steps, 6)];
				end
				%% Average probabilities for multiple steps in time boundaries
				steps = (0:(num_state_steps-1))*time_step/(num_state_steps-1);
				bound_up = temp2_table(:, 4) - this_time - steps;
				bound_lo = temp2_table(:, 3) - this_time - steps;
				distribution = vertex_clusters{this_vertex, next_vertexes(j), 3};
				try
				p_up = cdf(distribution{2}, bound_up - distribution{1}, distribution{3:end});
				catch
					bound_up
				end
				p_lo = cdf(distribution{2}, bound_lo - distribution{1}, distribution{3:end});
				prob = sum(p_up - p_lo, 2) / num_state_steps;
				temp_transition_table = [states(i, 1) * ones(size(temp2_table, 1), 1), next_vertexes(j) * ones(size(temp2_table, 1), 1), temp2_table(:, 1), prob, -1 * vertex_cluster_rewards(this_vertex) * ones(size(temp2_table, 1), 1), zeros(size(temp2_table, 1), 1)];
				idx = (temp_transition_table(:, 4) > 0);
				state_transition_table(counter + [1:sum(idx)], :) = temp_transition_table(idx, :);
				state_transition_table(counter + sum(idx) + 1, :) = [states(i ,1), next_vertexes(j), failure_states(this_vertex), 1 - sum(prob), -1 * vertex_cluster_rewards(this_vertex), 0];
				counter = counter + sum(idx) + 1;
			end
		elseif this_vertex == num_vertex
			state_transition_table(counter + 1, :) = [states(i, 1), absorbing_vertex, absorbing_state, 1, -1 * vertex_cluster_rewards(this_vertex), 0];
			counter = counter + 1;
		elseif this_vertex < absorbing_vertex
			state_transition_table(counter + 1, :) = [states(i, 1), absorbing_vertex, absorbing_state, 1, failure_penalties(states(i, 2) - num_vertex), 1];
			counter = counter + 1;
		else
			state_transition_table(counter + 1, :) = [states(i, 1), absorbing_vertex, absorbing_state, 1, 0, 0];
			counter = counter + 1;
		end
	end
	state_transition_table(counter+1:end, :) = [];
	%% Convert tables
	states = array2table(states, 'VariableNames', {'state', 'vertex', 'start_time', 'end_time'});
	state_transition_table = state_transition_table(state_transition_table(:, 4) > 0, :);
	state_transition_table = array2table(state_transition_table, 'VariableNames', {'state', 'maneuver', 'next_state', 'probability', 'value_0', 'value_1'});

end
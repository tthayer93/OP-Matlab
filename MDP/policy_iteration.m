function [values, policy, evals] = policy_iteration(state_transition_table, discount_rate, threshold, max_iterations, initial_policy)
%POLICY_ITERATION Compute an MDP policy using policy iteration
%
%	Version: 2.0
%	Date: 11/06/2020
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function uses the policy iteration algorithm to compute a complete policy for a Marko Decision Process (MDP) defined by the state transition table with the given discount rate. Every iteration computes state values and evals in-place, in the order that states appear in state_transition_table.
%	Inputs:
%		state_transition_table: a table containing all possible state transitions (one per row) and associated costs, with the following format
%			state | maneuver | next_state | probability | cost
%				All values in the state_transition_table must be numerals. If extra columns are provided, they will be used for additional evaluations of each state for the policy, without effecting the policy itself.
%		discount_rate: a value between zero and one defining how to weight future costs
%		threshold: a value close to zero defining the stopping critera for the algorithm
%		max_iterations: an integer defining the max number of iterations for the algorithm
%		initial_policy: a table containing the policy on which to begin iterative improvement, in the following format
%			state | maneuver
%	Outputs:
%		values: a table containing the computed values for each state, in the following format
%			state | value
%		policy: a table containing the computed policy for the MDP, in the following format
%			state | maneuver
%		evals: if extra columns in the state_transition_table are provided, this is a matrix containing extra evaluations for each state, where each row corresponds to the state provided in the same row within values.

	
	%% Initialize
	if ~exist('state_transition_table')
		error('Must provide state_transition_table.');
	else
		if size(state_transition_table, 2) < 5
			error('state_transition_table must have at least 5 columns.');
		end
	end
	if ~exist('discount_rate') || isempty(discount_rate)
		discount_rate = 1;
	else
		if discount_rate < 0 || discount_rate > 1
			error('discount_rate must be between 0 and 1.');
		end
	end
	if ~exist('threshold') || isempty(threshold)
		threshold = 0.001;
	else
		if threshold < 0
			error('threshold must be greater than 0.');
		end
	end
	if ~exist('max_iterations') || isempty(threshold)
		max_iterations = 1000;
	else
		if max_iterations < 0
			error('max_iterations must be greater than 0.');
		end
	end
	states = unique(state_transition_table(:, 1), 'stable');
	values = zeros(length(states), size(state_transition_table, 2) - 4);
	policy = [states, zeros(length(states), 1)];
	%% Initialize policy
	if ~exist('initial_policy') || isempty(initial_policy)
		for i=1:length(states)
			idxs = (state_transition_table(:, 1) == states(i));
			this_state_table = state_transition_table(idxs, :);
			temp_maneuvers = unique(this_state_table(:, 2));
			policy(i, :) = [states(i), temp_maneuvers(randperm(length(temp_maneuvers), 1))];
		end
	else
		if size(initial_policy, 1) ~= length(states) || size(initial_policy, 2) ~= 2
			error('initial_policy given in incorrect format.');
		else
			policy = initial_policy;
		end
	end	
	%% Start policy iteration
	policy_stable = false;
	while ~policy_stable
		%% Policy evaluation
		delta = inf;
		iteration = 0;
		while delta > threshold && iteration <= max_iterations
			iteration = iteration + 1;
			delta = 0;
			%% Loop over every state and get value for policy
			for i=1:length(states)
				prev_value = values(i, 1);
				this_maneuver = policy(policy(:, 1) == states(i), 2);
				idxs = (state_transition_table(:, 1) == states(i)) & (state_transition_table(:, 2) == this_maneuver);
				this_maneuver_table = state_transition_table(idxs, :);
				next_state_values = values(any(states(:, 1) == this_maneuver_table(:, 3)', 2), :);
				values(i, :) = sum(this_maneuver_table(:, 4) .* (this_maneuver_table(:, 5:end) + discount_rate * next_state_values), 1);
				delta = max(delta, abs(prev_value - values(i, 1)));
			end
		end
		%% Policy improvement
		policy_stable = true;
		for i=1:length(states)
			old_maneuver = policy(policy(:, 1) == states(i), 2);
			idxs = (state_transition_table(:, 1) == states(i));
			this_state_table = state_transition_table(idxs, :);
			temp_maneuvers = unique(this_state_table(:, 2));
			maneuver_values = zeros(length(temp_maneuvers), 1);
			%% Loop over maneuvers for current state
			for j=1:length(maneuver_values)
				this_maneuver_table = this_state_table(this_state_table(:, 2) == temp_maneuvers(j), :);
				if abs(1 - sum(this_maneuver_table(:, 4))) > 0.00001
					error('Transition probabilities for state %d and maneuver %d do not sum to 1.', states(i), temp_maneuvers(j));
				end
				next_state_values = values(any(states(:, 1) == this_maneuver_table(:, 3)', 2), 1);
				maneuver_values(j) = sum(this_maneuver_table(:, 4) .* (this_maneuver_table(:, 5) + discount_rate * next_state_values));
			end
			%% Get best maneuver for current state
			%[unused, temp] = max(maneuver_values(:, 1));
			[unused, temp] = min(maneuver_values(:, 1));
			policy(i, 2) = temp_maneuvers(temp);
			if old_maneuver ~= policy(i, 2)
				policy_stable = false;
			end
		end
	end
	evals = values(:, 2:end);
	values = [states, values(:, 1)];
end
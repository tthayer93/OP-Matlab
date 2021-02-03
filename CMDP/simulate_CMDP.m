function [state_path, maneuvers, values, failure] = simulate_CMDP(policy, state_action, state_transition_table, initial_state_distribution, absorbing_states, constraints)
%simulate_CMDP Simulates a single run of the policy for a CMDP
%
%	Version: 1.0
%	Date: 16/10/19
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function simulates a single run of the policy for a Constrained Markov Decision Process (CMDP).
%   Inputs:
%		policy: a vector containing the optimal policy, which follows the corrresponding state_action
%		state_action: A matrix containing all possible state-maneuver pairs, with the following format
%			state | maneuver
%		state_transition_table: a table containing all possible state transitions (one per row) and associated values where V >= 0, with the following format
%			state | maneuver | next_state | probability | value_0 | ... | value_V
%		initial_state_distribution: a matrix containing the probability of starting in each state, with the following format
%			state | probability
%		absorbing_states: a vector containing all absorbing states in the transition system
%		constraints: a vector containing the constraints of the CMDP, the length of which is less than or equal to V
%	Outputs:
%		state_path: a vector containing the sequence of states visited by the simulation
%		maneuvers: a vector containing the sequence of maneuvers chosen by the policy in simulation
%		values: a vector containing the sum of values for visiting each state, according to the state transition table
%		failure: a binary variable that is true if any of the constraints are violated

	%% Simulate the policy
	state_transition_table = table2array(state_transition_table);
	values = zeros(1, 1+length(constraints));
	maneuvers = zeros(length(policy), 1);
	maneuvers_counter = 0;
	state_path = zeros(length(policy), 1);
	state_path(1) = randsample(initial_state_distribution(:, 1), 1, 1, initial_state_distribution(:, 2));
	state_path_counter = 1;
	while ~any(state_path(state_path_counter) == absorbing_states)
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
		else
			maneuvers(maneuvers_counter + 1) = state_action(idx, 2);
			maneuvers_counter = maneuvers_counter + 1;
		end
		temp_table = state_transition_table((state_transition_table(:, 1) == state_path(state_path_counter)) & (state_transition_table(:, 2) == maneuvers(maneuvers_counter)), :);
		%% Get next state and update values
		if size(temp_table, 1) > 1
			state_path(state_path_counter + 1) = randsample(temp_table(:, 3), 1, 1, temp_table(:, 4));
			state_path_counter = state_path_counter + 1;
		else
			state_path(state_path_counter + 1) = temp_table(1, 3);
			state_path_counter = state_path_counter + 1;
		end
		values = values + temp_table(temp_table(:, 3) == state_path(state_path_counter), 5:end);
	end
	%% Shrink data structures to fit
	maneuvers((maneuvers_counter + 1):end) = [];
	state_path((state_path_counter + 1):end) = [];
	
	%% Check for constraint violation
	if exist('constraints')
		if any(values(2:end) > constraints)
			failure = 1;
		else
			failure = 0;
		end
	else
		failure = [];
	end

end


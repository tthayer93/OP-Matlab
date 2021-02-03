function [f, a, b, a_eq, b_eq, state_action] = CMDP_to_LP(state_transition_table, initial_state_distribution, absorbing_states, constraints)
%SOLVE_CMDP Formulates a generic CMDP as a LP
%
%	Version: 1.0
%	Date: 03/10/19
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function formulates a generic Constrained Markov Decision Process (CMDP) as a Linear Program (LP), such that the sum of rho * value_0 is minimized and the sums of rho * value_X are less than or equal to the constraint X.
%   Inputs:
%		state_transition_table: a table containing all possible state transitions (one per row) and associated values where N >= 0, with the following format
%			state | maneuver | next_state | probability | value_0 | ... | value_N
%		initial_state_distribution: a matrix containing the probability of starting in each state, with the following format
%			state | probability
%		absorbing_states: a vector containing all absorbing states in the transition system
%		constraints: a vector containing the constraints of the CMDP, the length of which is less than or equal to N
%	Outputs:
%		f: the function being minimized
%		a: the matrix of inequality functions
%		b: the matrix of inequality constraints
%		a_eq: the matrix of equality functions
%		b_eq: the matrix of equality constraints
%		state_action: A matrix containing all possible state-maneuver pairs, with the following format
%			state | maneuver
	
	epsilon = 0.00000001;
	%% Check state transition table and remove transitions with 0 probability
	if istable(state_transition_table)
		state_transition_table = table2array(state_transition_table);
	elseif ~ismatrix(state_transition_table)
		error('Unknown type for state_transition_table.');
	end
	state_transition_table(state_transition_table(:, 4) == 0, :) = [];
	%% Check initial state distribution
	if abs(1 - sum(initial_state_distribution(:, 2))) >= epsilon
		error('Sum of initial state probabilities must equal 1.');
	end
	%% Remove absorbing states
	states = unique([state_transition_table(:, 1); state_transition_table(:, 3)]);
	if exist('absorbing_states')
		absorbing_states = absorbing_states(:);
		states(any(states == absorbing_states', 2)) = [];
	end
	%% Check constraints
	if exist('constraints')
		b = constraints(:);
		if (size(state_transition_table, 2) - 5) < length(constraints)
			error('Too many constraints.');
		end
	else
		b = [];
	end
	%% Build linprog inputs
	state_action = zeros(size(state_transition_table, 1), 2);
	state_action_counter = 0;
	f = zeros(size(state_transition_table, 1), 1);
	f_counter = 0;
	a_list = zeros(size(state_transition_table, 1)*length(b), 3);
	a_counter = 0;
	aeq_list = zeros(size(state_transition_table, 1), 3);
	aeq_counter = 0;
	b_eq = zeros(size(state_transition_table, 1), 1);
	beq_counter = 0;
	for i=1:length(states)
		temp_table = state_transition_table(state_transition_table(:, 1) == states(i), :);
		maneuvers = unique(temp_table(:, 2));
		for j=1:length(maneuvers)
			state_action(state_action_counter + 1, :) = [states(i), maneuvers(j)];
			state_action_counter = state_action_counter + 1;
			temp2_table = temp_table(temp_table(:, 2) == maneuvers(j), :);
			if any(temp2_table(:, 5:end) ~= temp2_table(1, 5:end))
				error('All values for a state-maneuver pair in state_transition_table must be the same for each outcome.');
			end
			if abs(1 - sum(temp2_table(:, 4))) >= epsilon
				error('Sum of outcome probabilities for state-maneuver pair must equal 1.');
			end
			f(f_counter + 1) = temp2_table(1, 5);
			f_counter = f_counter + 1;
			this_a = temp2_table(1, 6:(5+length(b)))';
			nonzeros_a = find(this_a);
			a_list(a_counter + [1:length(nonzeros_a)], :) = [nonzeros_a, state_action_counter * ones(length(nonzeros_a), 1), this_a(nonzeros_a)];
			a_counter = a_counter + length(nonzeros_a);
			absorbing_idx = false(size(temp2_table, 1), 1);
			for k=1:length(absorbing_states)
				absorbing_idx = absorbing_idx | (temp2_table(:, 3) == absorbing_states(k));
			end
			temp3_table = temp2_table(~absorbing_idx, :);
			idxs = (temp3_table(:, 3)' == states);
			idxs1 = any(idxs, 1);
			idxs2 = any(idxs, 2);
			aeq_list(aeq_counter + 1, :) = [i, state_action_counter, 1];
			aeq_counter = aeq_counter + 1;
			aeq_list(aeq_counter + [1:sum(idxs1)], :) = [find(idxs2), state_action_counter * ones(sum(idxs1), 1), 0 - temp3_table(idxs1, 4)];
			aeq_counter = aeq_counter + sum(idxs1);
			if aeq_counter == size(aeq_list, 1)
				aeq_list = [aeq_list; zeros(size(state_transition_table, 1), 3)];
			end
		end
		idx = find(initial_state_distribution(:, 1) == states(i));
		if isempty(idx)
			b_eq(beq_counter + 1) = 0;
			beq_counter = beq_counter + 1;
		else
			b_eq(beq_counter + 1) = initial_state_distribution(idx, 2);
			beq_counter = beq_counter + 1;
		end
	end
	state_action((state_action_counter + 1):end, :) = [];
	f((f_counter + 1):end) = [];
	a_list((a_counter + 1):end, :) = [];
	a = sparse(a_list(:, 1), a_list(:, 2), a_list(: ,3), length(b), size(state_action, 1));
	b_eq((beq_counter + 1):end) = [];
	aeq_list((aeq_counter + 1):end, :) = [];
	a_eq = sparse(aeq_list(:, 1), aeq_list(:, 2), aeq_list(:, 3), length(b_eq), size(state_action, 1));

end
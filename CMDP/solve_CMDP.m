function [policy, state_action, rho, values] = solve_CMDP(state_transition_table, initial_state_distribution, absorbing_states, constraints, linprog_options)
%SOLVE_CMDP Solves a generic CMDP
%
%	Version: 2.01
%	Date: 11/06/20
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves a generic Constrained Markov Decision Process (CMDP), such that the sum of rho * value_0 is minimized and the sums of rho * value_X are less than or equal to the constraint X.
%   Inputs:
%		state_transition_table: a table containing all possible state transitions (one per row) and associated values where N >= 0, with the following format
%			state | maneuver | next_state | probability | value_0 | ... | value_N
%		initial_state_distribution: a matrix containing the probability of starting in each state, with the following format
%			state | probability
%		absorbing_states: a vector containing all absorbing states in the transition system
%		constraints: a vector containing the constraints of the CMDP, the length of which is less than or equal to N
%		linprog_options: pass through options for Matlab's linprog(), which if exists, will force usage of matlabs linear program solver
%	Outputs:
%		policy: a vector containing the optimal policy, which follows the corrresponding state_action
%		state_action: A matrix containing all possible state-maneuver pairs, with the following format
%			state | maneuver
%		rho: a vector containing the computed rho values of the CMDP
%		values: a vector containing the expected values of the CMDP
	
	%% Create LP from CMDP
	[f, a, b, a_eq, b_eq, state_action] = CMDP_to_LP(state_transition_table, initial_state_distribution, absorbing_states, constraints);
	%% Solve linear program
	lb = zeros(numel(f),1);
	ub = Inf*ones(numel(f),1);
	if exist('linprog_options')
		rho = solve_LP(f, a, b, a_eq, b_eq, lb, ub, linprog_options);
	else
		rho = solve_LP(f, a, b, a_eq, b_eq, lb, ub);
	end
	%% Compute policy and values
	[policy, values] = compute_CMDP_policy(rho, state_action, f, a);

end
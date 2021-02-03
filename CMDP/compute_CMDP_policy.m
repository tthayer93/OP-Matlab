function [policy, values] = compute_CMDP_policy(rho, state_action, f, a)
%SOLVE_CMDP Computes the policy of a CMDP
%
%	Version: 1.1
%	Date: 03/10/19
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function computes the policy of a Constrained Markov Decision Process (CMDP) from a given rho vector, which is the output of a linear program solving the CMDP.
%   Inputs:
%		rho: a vector containing the computed rho values of the CMDP
%		state_action: A matrix containing all possible state-maneuver pairs, with the following format
%			state | maneuver
%		a: the matrix of inequality functions for the linear program solving the CMDP
%	Outputs:
%		policy: a vector containing the optimal policy, which follows the corrresponding state_action
%		values: a vector containing the expected values of the CMDP
	
	%epsilon = 0.00000001;
	epsilon = 0;
	%% Compute policy
	temp_state = state_action(1, 1);
	temp_i = 1;
	summ = 0;
	policy = zeros(length(rho), 1);
	for i=1:size(state_action, 1)
		if (state_action(i, 1) ~= temp_state) & (i > 1)
			policy((i-1):-1:temp_i) = policy((i-1):-1:temp_i) / summ;
			temp_state = state_action(i, 1);
			temp_i = i;
			summ = 0;
		end
		if rho(i) > epsilon
			policy(i) = rho(i);
			summ = summ + rho(i);
		else
			policy(i) = 0;
		end
	end
	policy(isnan(policy)) = 0;
	policy(policy == inf) = 1;
	%% Compute values
	if exist('f') & exist('a')
		values = zeros(size(a, 1)+1, 1);
		values(1) = dot(rho, f);
		for row_index = 1:size(a, 1)
			values(row_index + 1) = dot(rho,a(row_index, :));
		end
	else
		values = [];
	end

end
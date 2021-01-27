function [tour, total_reward, total_cost] = solve_OP(edge_list, vertex_rewards, budget, beginning, ending, intlinprog_options)
%SOLVE_OP Solves a generic OP
%
%	Version: 1.0
%	Date: 25/11/19
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves a generic Orienteering Problem (OP), such that the sum of rewards is maximized for the given budget constraint.
%	Inputs:
%		edge_list: a list containing the set of directed edges in the graph, with the following format
%			initial vertex | destination vertex | distance
%		vertex_rewards: a vector containing the reward for each vertex
%		budget: a scalar bound on the total distance of the computed route
%		beginning: the first vertex in the route
%		ending: the last vertex in the route
%		intlinprog_options: pass through options for Matlab's linprog(), which if exists, will force usage of matlabs linear program solver
%	Outputs:
%		tour: a vector containing the sequence of verticies visited by the optimal tour
%		total_reward: the total reward collected by the computed tour
%		total_cost: the total cost associated with the computed tour

	%% Create ILP from OP
	[f, x_int, a, b, a_eq, b_eq, lb, ub] = OP_to_ILP(edge_list, vertex_rewards, budget, beginning, ending);
	%% Solve ILP
	if exist('intlinprog_options')
		[x_opt] = solve_ILP(f, x_int, a, b, a_eq, b_eq, lb, ub, intlinprog_options);
	else
		[x_opt] = solve_ILP(f, x_int, a, b, a_eq, b_eq, lb, ub);
	end
	x_opt = x_opt(1:size(edge_list, 1));
	%% Compute tour
	tour = [beginning];
	total_reward = vertex_rewards(beginning);
	vertex_rewards(beginning) = 0;
	total_cost = 0;
	while 1
		edge = find(round(x_opt) & edge_list(:, 1) == tour(end));
		tour(end+1) = edge_list(edge, 2);
		total_reward = total_reward + vertex_rewards(tour(end));
		vertex_rewards(tour(end)) = 0;
		total_cost = total_cost + edge_list(edge, 3);
		if tour(end) == ending
			break;
		end
	end
	
end


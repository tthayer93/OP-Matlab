function [tour, total_cost] = solve_TSP(edge_list, symmetric, intlinprog_options)
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
%		symmetric: a boolean variable describing if edges are symmetric
%		linprog_options: pass through options for Matlab's linprog(), which if exists, will force usage of matlabs linear program solver
%	Outputs:
%		tour: a vector containing the sequence of verticies visited by the optimal tour
%		total_cost: the total cost associated with the computed tour

	%% Initialize
	if ~exist('symmetric') || isempty(symmetric)
		symmetric = 0;
	end
	%% Create ILP from TSP
	[f, x_int, a, b, a_eq, b_eq, lb, ub] = TSP_to_ILP(edge_list, symmetric);
	%% Solve ILP
	if exist('intlinprog_options')
		[x_opt] = solve_ILP(f, x_int, a, b, a_eq, b_eq, lb, ub, intlinprog_options);
	else
		[x_opt] = solve_ILP(f, x_int, a, b, a_eq, b_eq, lb, ub);
	end
	x_opt = x_opt(1:size(edge_list, 1));
	%% Compute tour
	tour = [edge_list(1, 1)];
	total_cost = 0;
	while 1
		if symmetric
			edge = find(round(x_opt) & (edge_list(:, 1) == tour(end) | edge_list(:, 2) == tour(end)));
			if length(edge) > 1
				edge = edge(1);
			end
			if edge_list(edge, 1) == tour(end)
				tour(end+1) = edge_list(edge, 2);
			else
				tour(end+1) = edge_list(edge, 1);
			end
			total_cost = total_cost + edge_list(edge, 3);
			edge_list(edge, :) = [];
			x_opt(edge, :) = [];
		else
			edge = find(round(x_opt) & edge_list(:, 1) == tour(end));
			tour(end+1) = edge_list(edge, 2);
			total_cost = total_cost + edge_list(edge, 3);
		end
		if tour(end) == tour(1)
			break;
		end
	end
	
end


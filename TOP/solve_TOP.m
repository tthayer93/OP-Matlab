function [tours, total_rewards, total_costs] = solve_TOP(edge_list, vertex_rewards, budgets, beginning, ending, intlinprog_options)
%SOLVE_TOP Solves a generic TOP
%
%	Version: 1.0
%	Date: 25/11/19
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves a generic Team Orienteering Problem (TOP), such that the sum of rewards is maximized for the given budget constraints of each agent.
%	Assumptions:
%		The graph should be fully connected, or tour extraction will not work for all agents.
%	Inputs:
%		edge_list: a list containing the set of directed edges in the graph, with the following format
%			initial vertex | destination vertex | distance
%		vertex_rewards: a vector containing the reward for each vertex
%		budgets: a vector containing the bounds on the total distance of the computed route for each agent
%		beginning: the first vertex in each route
%		ending: the last vertex in each route
%		intlinprog_options: pass through options for Matlab's linprog(), which if exists, will force usage of matlabs linear program solver
%	Outputs:
%		tour: a vector containing the sequence of verticies visited by the optimal tour
%		total_reward: the total reward collected by the computed tour
%		total_cost: the total cost associated with the computed tour

	%% Create ILP from OP
	[f, x_int, a, b, a_eq, b_eq, lb, ub] = TOP_to_ILP(edge_list, vertex_rewards, budgets, beginning, ending);
	%% Solve ILP
	if exist('intlinprog_options')
		[x_opt] = solve_ILP(f, x_int, a, b, a_eq, b_eq, lb, ub, intlinprog_options);
	else
		[x_opt] = solve_ILP(f, x_int, a, b, a_eq, b_eq, lb, ub);
	end
	x_opt = x_opt(1:size(edge_list, 1));
	%% Compute tours
	tours = cell(1, length(budgets));
	for i = 1:length(budgets)
		this_tour = [beginning];
		total_rewards = vertex_rewards(beginning);
		vertex_rewards(beginning) = 0;
		total_costs = 0;
		while 1
			edge = find(round(x_opt) & edge_list(:, 1) == this_tour(end), 1);
			x_opt(edge) = 0;
			this_tour(end+1) = edge_list(edge, 2);
			total_rewards = total_rewards + vertex_rewards(this_tour(end));
			vertex_rewards(this_tour(end)) = 0;
			total_costs = total_costs + edge_list(edge, 3);
			if this_tour(end) == ending
				break;
			end
		end
		tours{i} = this_tour;
	end
	
end
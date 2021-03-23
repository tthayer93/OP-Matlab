function [f, x_int, a, b, a_eq, b_eq, lb, ub] = TOP_to_ILP(edge_list, vertex_rewards, budgets, beginning, ending)
%TOP_TO_ILP Formulates a generic TOP as an ILP
%
%	Version: 1.0
%	Date: 23/03/21
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function formulates a generic Team Orienteering Problem (OP) as an Integer Linear Program (ILP), such that the sum of f * x is minimized.
%	Inputs:
%		edge_list: a list containing the set of directed edges in the graph, with the following format
%			initial vertex | destination vertex | distance
%		vertex_rewards: a vector containing the reward for each vertex
%		budgets: a vector containing the bounds on the total distance of the computed route for each agent
%		beginning: the first vertex in the route
%		ending: the last vertex in the route
%	Outputs:
%		f: the function being minimized
%		x_int: the x variables that are contrained as integers
%		a: the matrix of inequality functions
%		b: the vector of inequality constraints
%		a_eq: the matrix of equality functions
%		b_eq: the vector of equality constraints
%		lb: a vector containing the lower bound for each x
%		ub: a vector containing the upper bound for each x
%
% TODO:
%	CURRENTLY NOT WORKING! DO NOT USE! No solution satisfies constraints
%	Different cost for each agent
%	Different reward for each agent
%	Different start and end vertices for each agent

	%% Setup f
	n_vertex = max([length(vertex_rewards); edge_list(:, 1); edge_list(:, 2)]);
	n_agents = length(budgets);
	rewards = zeros(1, n_vertex);
	rewards(1:length(vertex_rewards)) = vertex_rewards;
	agent_edge_list = repmat(edge_list, n_agents, 1);
	f = [-rewards(agent_edge_list(:, 2)), zeros(1, n_vertex*n_agents)];
	%% Setup inequality matrix and constraints
	a = [];
	for x=1:n_agents
		these_edges = [zeros((x-1)*size(edge_list, 1), 1); ones(size(edge_list, 1), 1); zeros((n_agents-x)*size(edge_list, 1), 1)];
		a = [a; these_edges'.*agent_edge_list(:, 3)'.*spones(1:size(agent_edge_list, 1))];
	end
	b = budgets(:);
	a = [a; sparse(n_vertex, size(agent_edge_list, 1))];
	k = 1;
	for i=1:n_vertex
		k = k + 1;
		leaving_indexes = (agent_edge_list(:, 1) == i);
		arriving_indexes = (agent_edge_list(:, 2) == i);
		indexes = sparse(sum(leaving_indexes + arriving_indexes, 2));
		a(k, :) = indexes';
		if (i == beginning) || (i == ending)
			b = [b; 2*n_agents];
		else
			b = [b; 2];
		end
	end
	%% Setup Miller Tucker Zemlin constraints
	a = [a, sparse(size(a, 1), n_vertex*n_agents)];
	temp = repmat(1:n_vertex, 1, n_agents);
	for x=1:n_agents
		these_edges = [zeros((x-1)*size(edge_list, 1), 1); ones(size(edge_list, 1), 1); zeros((n_agents-x)*size(edge_list, 1), 1)];
		for i=1:n_vertex
			for j=2:n_vertex
				mtz = (these_edges .* (agent_edge_list(:, 1) == i) .* (agent_edge_list(:, 2) == j) * (n_vertex - 1))';
				mtz = [mtz, (temp == i) - (temp == j)];
				a = [a; sparse(mtz)];
				b = [b; n_vertex - 2];
			end
		end
	end
	%% Setup equality constraints
	if beginning ~= ending
		a_eq = 1.0*[((agent_edge_list(:, 1) == beginning) - (agent_edge_list(:, 2) == beginning))'; ((agent_edge_list(:, 1) == ending) - (agent_edge_list(:, 2) == ending))'];
		b_eq = [n_agents; -n_agents];
		a_eq = [sparse(a_eq); sparse(n_vertex - 2, size(agent_edge_list, 1))];
		k = 2;
	else
		a_eq = 1.0*[((agent_edge_list(:, 1) == beginning) - (agent_edge_list(:, 2) == beginning))'; ((agent_edge_list(:, 1) == beginning) + (agent_edge_list(:, 2) == beginning))'];
		b_eq = [0; 2*n_agents];
		a_eq = [sparse(a_eq); sparse(n_vertex - 1, length(agent_edge_list))];
		k = 2;
	end
	for x=1:n_agents
		these_edges = [zeros((x-1)*size(edge_list, 1), 1); ones(size(edge_list, 1), 1); zeros((n_agents-x)*size(edge_list, 1), 1)];
		for i=1:n_vertex
			if i~=beginning && i~=ending
				k = k + 1;
				a_eq(k, :) = (these_edges.*(agent_edge_list(:, 1) == i) - (these_edges.*(agent_edge_list(:, 2) == i)))';
				b_eq = [b_eq; 0];
			end
		end
	end
	a_eq = [a_eq, sparse(size(a_eq, 1), n_vertex*n_agents)];
	a_eq = [a_eq; [zeros(1, size(agent_edge_list, 1)), temp == beginning]];
	b_eq = [b_eq; 1];
	%% Setup bounds
	x_int = 1:(size(agent_edge_list, 1) + n_vertex*n_agents);
	lb = zeros(1, (size(agent_edge_list, 1) + n_vertex*n_agents));
	ub = ones(1, size(agent_edge_list, 1));
	ub = [ub, n_vertex * ones(1, n_vertex*n_agents)];
end
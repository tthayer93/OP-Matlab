function [f, x_int, a, b, a_eq, b_eq, lb, ub] = OP_to_ILP(edge_list, vertex_rewards, budget, beginning, ending)
%OP_TO_ILP Formulates a generic OP as an ILP
%
%	Version: 1.0
%	Date: 25/11/19
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function formulates a generic Orienteering Problem (OP) as an Integer Linear Program (ILP), such that the sum of f * x is minimized.
%	Inputs:
%		edge_list: a list containing the set of directed edges in the graph, with the following format
%			initial vertex | destination vertex | distance
%		vertex_rewards: a vector containing the reward for each vertex
%		budget: a scalar bound on the total distance of the computed route
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

	%% Setup f
	n_vertex = max([length(vertex_rewards); edge_list(:, 1); edge_list(:, 2)]);
	rewards = zeros(1, n_vertex);
	rewards(1:length(vertex_rewards)) = vertex_rewards;
	f = [-rewards(edge_list(:, 2)), zeros(1, n_vertex)];
	%% Setup inequality matrix and constraints
	a = edge_list(:, 3)'.*spones(1:size(edge_list, 1));
	b = budget;
	a = [a; sparse(n_vertex, size(edge_list, 1))];
	k = 1;
	for i=1:n_vertex
		k = k + 1;
		leaving_indexes = (edge_list(:, 1) == i);
		arriving_indexes = (edge_list(:, 2) == i);
		indexes = sparse(sum(leaving_indexes + arriving_indexes, 2));
		a(k, :) = indexes';
		b = [b; 2];
	end
	%% Setup Miller Tucker Zemlin constraints
	a = [a, sparse(size(a, 1), n_vertex)];
	temp = 1:n_vertex;
	for i=1:n_vertex
		for j=2:n_vertex
			mtz = ((edge_list(:, 1) == i) .* (edge_list(:, 2) == j) * (n_vertex - 1))';
			mtz = [mtz, (temp == i) - (temp == j)];
			a = [a; sparse(mtz)];
			b = [b; n_vertex - 2];
		end
	end
	%% Setup equality constraints
	if beginning ~= ending
		a_eq = 1.0*[((edge_list(:, 1) == beginning) - (edge_list(:, 2) == beginning))'; ((edge_list(:, 1) == ending) - (edge_list(:, 2) == ending))'];
		b_eq = [1; -1];
		a_eq = [sparse(a_eq); sparse(n_vertex - 2, size(edge_list, 1))];
		k = 2;
	else
		a_eq = 1.0*[((edge_list(:, 1) == beginning) - (edge_list(:, 2) == beginning))'; ((edge_list(:, 1) == beginning) + (edge_list(:, 2) == beginning))'];
		b_eq = [0; 2];
		a_eq = [sparse(a_eq); sparse(n_vertex - 1, length(edge_list))];
		k = 2;
	end
	for i=1:n_vertex
		if i~=beginning && i~=ending
			k = k + 1;
			a_eq(k, :) = (1.0*(edge_list(:, 1) == i) - (1.0*(edge_list(:, 2) == i)))';
			b_eq = [b_eq; 0];
		end
	end
	a_eq = [a_eq, sparse(size(a_eq, 1), n_vertex)];
	a_eq = [a_eq; [zeros(1, size(edge_list, 1)), temp == beginning]];
	b_eq = [b_eq; 1];
	%% Setup bounds
	x_int = 1:(size(edge_list, 1) + n_vertex);
	lb = zeros(1, (size(edge_list, 1) + n_vertex));
	ub = ones(1, size(edge_list, 1));
	ub = [ub, n_vertex * ones(1, n_vertex)];
	
end


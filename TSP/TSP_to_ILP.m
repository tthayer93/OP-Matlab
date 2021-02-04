function [f, x_int, a, b, a_eq, b_eq, lb, ub] = TSP_to_ILP(edge_list, symmetric)
%TSP_TO_ILP Formulates a TSP as an ILP
%
%	Version: 1.0
%	Date: 15/07/20
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function formulates a Traveling Salesman Problem (TSP) as an Integer Linear Program (ILP), such that the sum of f * x is minimized.
%	Inputs:
%		edge_list: a list containing the set of directed edges in the graph, with the following format
%			initial vertex | destination vertex | distance
%		symmetric: a boolean variable describing if edges are symmetric
%	Outputs:
%		f: the function being minimized
%		x_int: the x variables that are contrained as integers
%		a: the matrix of inequality functions
%		b: the vector of inequality constraints
%		a_eq: the matrix of equality functions
%		b_eq: the vector of equality constraints
%		lb: a vector containing the lower bound for each x
%		ub: a vector containing the upper bound for each x
	
	%% Initialize
	if ~exist('symmetric') || isempty(symmetric)
		symmetric = 0;
	end
	if symmetric
		%test if actually symmetric
	end
	%% Setup f
	vertices = unique([edge_list(:, 1); edge_list(:, 2)]);
	n_vertex = length(vertices);
	f = edge_list(:, 3)';
	%% Setup equality constraints
	a_eq = sparse(0, size(edge_list, 1));
	b_eq = [];
	for i=1:n_vertex
		if symmetric
			a_eq = [a_eq; (1.0*(edge_list(:, 1) == vertices(i) | edge_list(:, 2) == vertices(i)))'];
			b_eq = [b_eq; 2];
		else
			a_eq = [a_eq; (1.0*(edge_list(:, 1) == vertices(i)))'];
			b_eq = [b_eq; 1];
			a_eq = [a_eq; (1.0*(edge_list(:, 2) == vertices(i)))'];
			b_eq = [b_eq; 1];
		end
	end
	%% Setup inequality constraints
	a = sparse(0, size(edge_list, 1));
	b = [];
% 	for i = 1:n_vertex
% 		for j = 1:n_vertex
% 			for s=2:(n_vertex-2)
% 				a = [a; (edge_list(:, 1) == vertices(i) & edge_list(:, 2) == verticies(j))'];
% 				b = [b; s-1];
% 			end
% 		end
% 	end
	%% Setup bounds
	x_int = 1:(size(edge_list, 1));
	lb = zeros(1, (size(edge_list, 1)));
	ub = ones(1, size(edge_list, 1));
	%% Setup Miller-Tucker-Zemlin constraints
	f = [f, zeros(1, n_vertex)];
	a_eq = [a_eq, sparse(size(a_eq, 1), n_vertex)];
	a = [a, sparse(size(a, 1), n_vertex)];
	x_int = [x_int, max(x_int+1):max(x_int+n_vertex)];
	lb = [lb, zeros(1, n_vertex)];
	ub = [ub, (n_vertex - 1) * ones(1, n_vertex)];
	for i=1:n_vertex
		for j=2:n_vertex
			if i ~= j
				mtz = ((edge_list(:, 1) == vertices(i)) .* (edge_list(:, 2) == vertices(j)) * n_vertex)';
				mtz = [mtz, (vertices' == vertices(i)) - (vertices' == vertices(j))];
				a = [a; sparse(mtz)];
				b = [b; n_vertex - 1];
			end
		end
	end
	
end


function [vertex_cluster_rewards, vertex_clusters] = OP_route_to_SCOP(tour, edge_list, rewards, cost_dist_type, cost_min_frac)
%OP_ROUTE_TO_SCOP Creates a SCOP from an OP route
%
%	Version: 1.0
%	Date: 12/02/20
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function takes a route built to solve the Orienteering Problem (OP) creates a clustered Stochastic Cost Orienteering Problem (SCOP) from the route, utilizing the given edge cost as the expected cost. First used in https://ieeexplore.ieee.org/abstract/document/9340899
%	Assumptions:
%		Vertices are not clustered, but the data structure used to store SCOP information is the same
%	Inputs:
%		tour: the initial OP tour to base the SCOP on
%		edge_list: the available directed edges for travel in the OP, with the following format
%			starting_vertex | ending_vertex | expected_travel_cost
%		rewards: a list containing the rewards for visiting vertices, with the following format
%			vertex | reward
%		cost_dist_type: the type of single parameter distribution used to model costs between adjacent vertices in the vineyard
%		cost_min_frac: the fraction of cost for the minimum bound of each possible transition, such that cost*cost_min_frac + distribution(1-cost_min_frac) = cost of transition
%	Outputs:
%		vertex_cluster_rewards: the reward associated with visiting each vertex cluster
%		vertex_clusters: an upper triangular cell matrix of size NxNx3 that contains the following information
%			cell (x,y,1): a vector of the sequence of verticies to traverse from the end of the current cluster x (row) to the end of cluster y (column)
%			cell (x,y,2): a matrix (whos lengths match the corresponding sequence of verticies) of cells each containing
%				1: the minimum cost for the associated vertex transition
%				2: the type of single parameter distribution used to model the cost for the associated vertex transition
%				3: the parameter of the distribution for the associated vertex transition
%			cell (x,y,3): a cell containing the combined distribution of the sequence of verticies for the associated vertex cluster transition with
%				1: the minimum cost for the associated vertex cluster transition
%				2: the type of distribution used to model the cost for the associated vertex cluster transition
%				3+: the parameters of the distribution for the associated vertex cluster transition
%
%TODO
% shortest path between vertices in tour
% more than just gamma

	%% Check inputs
	if iscell(cost_dist_type)
		if size(cost_dist_type, 1) == 1 && size(cost_dist_type, 2) == 1
			temp = cost_dist_type;
			cost_dist_type = cell(size(edge_list, 1), 1);
			cost_dist_type(:) = {temp};
		end
		if max(size(cost_dist_type)) ~= size(edge_list, 1) || min(size(cost_dist_type)) ~= 1
			error('Length of cost_dist_type must equal the number of edges in edge_list.');	
		end
	else
		temp = cost_dist_type;
		cost_dist_type = cell(size(edge_list, 1), 1);
		cost_dist_type(:) = {temp};
	end
	for i=1:length(cost_dist_type)
		if ~ischar(cost_dist_type{i})
			error('cost_dist_type can contain only strings.');
		else
			if ~strcmp(cost_dist_type{i}, 'Exponential') && ~strcmp(cost_dist_type{i}, 'Gamma')
				error('Distribution type not supported.');
			end
		end
	end
	if max(size(cost_min_frac)) ~= size(edge_list, 1) || min(size(cost_min_frac)) ~= 1
		if max(size(cost_min_frac)) ~= 1
			error('cost_min_frac must be of length 1 or the same length as edge_list.');
		else
			cost_min_frac = cost_min_frac * ones(size(edge_list, 1), 1);
		end
	end
	%% Build vertex clusters
	%vertex_cluster_rewards = rewards(any(rewards(:, 1) == tour, 2), 2);
	vertex_cluster_rewards = zeros(length(tour), 1);
	vertex_clusters = cell(length(tour), length(tour), 3);
	for i=1:length(tour)
		vertex_cluster_rewards(i) = rewards(rewards(:, 1) == tour(i), 2);
		for j=(i+1):length(tour)
			if tour(i) ~= tour(j)
				edge_idx = (edge_list(:, 1) == tour(i)) & (edge_list(:, 2) == tour(j));
				vertex_clusters{i, j, 1} = [tour(j)];
				vertex_clusters{i, j, 2} = [num2cell(cost_min_frac(edge_idx) * edge_list(edge_idx, 3)), cost_dist_type(edge_idx), num2cell(1 - cost_min_frac(edge_idx))];
				vertex_clusters{i, j, 3} = {cost_min_frac(edge_idx) * edge_list(edge_idx, 3), 'Gamma', edge_list(edge_idx, 3), 1 - cost_min_frac(edge_idx)};
			else
				vertex_clusters{i, j, 1} = [tour(j)];
				vertex_clusters{i, j, 2} = [num2cell(0), {'None'}, 0];
				vertex_clusters{i, j, 3} = {0, 'Gamma', 0, 0};
			end
		end
	end
	
end


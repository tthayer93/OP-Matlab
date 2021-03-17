function [vertex_cluster_rewards, vertex_clusters] = vineyard_route_to_SCOP(vine_distance, row_distance, reward_map, tour, cost_dist_type, cost_min_frac)
%VINEYARD_ROUTE_TO_SCOP Creates a SCOP from a vineyard route
%
%	Version: 1.0
%	Date: 29/10/19
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function takes a route built over a vineyard block and creates a clustered Stochastic Cost Orienteering Problem (SCOP) from the route, where each row or partial row is a cluster. First used in https://ieeexplore.ieee.org/abstract/document/9340899
%	Assumptions:
%		
%	Inputs:
%		vine_distance: the distance between each vine in a row
%		row_distance: the distance between each row in the vineyard block
%		reward_map: a matrix of size LxW containing the reward for visiting each vertex in the vineyard, where L is the number of rows and W is the number of vines per row
%		tour: the initial tour to base the SCOP on
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

	%% Extract basic information
	size_row = size(reward_map, 1);
	size_column = size(reward_map, 2);
	vertex_list = 1:size_row*size_column;
	left_side_verticies = vertex_list(mod(vertex_list, size_column) == 1);
	right_side_verticies = vertex_list(mod(vertex_list, size_column) == 0);
	%% Figure out which side we are on, left=1 and right=2
	which_side = zeros(length(tour), 1);
	for i=1:length(tour)
		if any(tour(i) == left_side_verticies)
			which_side(i) = 1;
		end
		if any(tour(i) == right_side_verticies)
			which_side(i) = 2;
		end
	end
	%% Logic
	i = 1;
	v = 1;
	vertex_sequences = {[tour(1)]};
	vertex_cluster_rewards = [0];
	expected_costs = [0];
	if which_side(i) == 1 % first vertex is considered a partial row
		left_partial_vertexes = [v];
		right_partial_vertexes = [];
	else
		left_partial_vertexes = [];
		right_partial_vertexes = [v];
	end
	left_full_vertexes = [];
	right_full_vertexes = [];
	while i <= length(tour)
		this_reward = 0;
		this_cost = 0;
		if i < length(tour) % if not at end of tour
			while which_side(i+1) ~= 0 % go up or down side columns
				i = i + 1;
				if i == length(tour) % if at ending position
					break;
				end
			end
			if i < length(tour)
				this_sequence = [tour(i)];
				this_cost = this_cost + vine_distance;
				j = 1;
				while which_side(i+j) == 0 % go into row
					this_sequence = [this_sequence; tour(i+j)];
					if ~any(this_sequence(end) == this_sequence(1:(end-1)))
						this_reward = this_reward + reward_map(floor((tour(i+j) - 1) / size_column) + 1, mod(tour(i+j) - 1, size_column) + 1);
						reward_map(floor((tour(i+j) - 1) / size_column) + 1, mod(tour(i+j) - 1, size_column) + 1) = 0;
					end
					this_cost = this_cost + vine_distance;
					j = j + 1;
				end
				this_sequence = [this_sequence; tour(i+j)];
			end
		end
		if i == length(tour) % if at end of tour
			this_sequence = [tour(i)];
			j = 0;
		end
		v = v + 1;
		vertex_sequences{v} = this_sequence;
		vertex_cluster_rewards(v) = this_reward;
		if which_side(i) == which_side(i+j) % if partial row
			if which_side(i) == 1 % if left partial
				selectables = [left_partial_vertexes; right_full_vertexes];
				left_partial_vertexes = [left_partial_vertexes; v];
			else % if right partial
				selectables = [right_partial_vertexes; left_full_vertexes];
				right_partial_vertexes = [right_partial_vertexes; v];
			end
		else % if full row
			if which_side(i) == 1 % if left full
				selectables = [right_full_vertexes; left_partial_vertexes];
				left_full_vertexes = [left_full_vertexes; v];
			else % if right full
				selectables = [left_full_vertexes; right_partial_vertexes];
				right_full_vertexes = [right_full_vertexes; v];
			end
		end
		costs_vec = zeros(v, 1);
		for k=1:v-1
			if any(k == selectables) % selectable is any vertex ending on the side this one begins on
				k_row = floor((vertex_sequences{k}(end) - 1) / size_column) + 1;
				v_row = floor((vertex_sequences{v}(1) - 1) / size_column) + 1;
				side_cost = row_distance*abs(v_row - k_row);
				costs_vec(k) = side_cost + this_cost;
			else
				costs_vec(k) = inf;
			end
		end
		expected_costs = [[expected_costs; zeros(1,v-1)], costs_vec];
		if j == 0
			j = 1;
		end
		i = i + j;
	end
	%% Build vertex clusters
	vertex_clusters = cell(v, v, 3);
	for i=1:length(vertex_sequences)
		for j=i+1:length(vertex_sequences)
			%% Build sequence of verticies
			if any(vertex_sequences{i}(end) == left_side_verticies) == any(vertex_sequences{j}(1) == left_side_verticies) % check end of i is on same side as beginning of j
				i_row = floor((vertex_sequences{i}(end) - 1) / size_column) + 1;
				j_row = floor((vertex_sequences{j}(1) - 1) / size_column) + 1;
				this_sequence = [];
				derp = 1;
				if i_row > j_row
					derp = -1;
				end
				for k=(i_row+derp):derp:j_row
					if any(vertex_sequences{i}(end) == left_side_verticies)
						this_sequence = [this_sequence; (k-1)*size_column+1];
					else
						this_sequence = [this_sequence; k*size_column];
					end
				end
				vertex_clusters{i, j, 1} = [this_sequence; vertex_sequences{j}(2:end)];
			end
			%% Build sequence of distributions
			this_sequence = vertex_clusters{i, j, 1};
			size_sequence = length(this_sequence);
			expected_costs = zeros(size_sequence, 1);
			if size_sequence >= 1
				if any(vertex_sequences{i}(end) == left_side_verticies) && any(this_sequence(1) == left_side_verticies)
					expected_costs(1) = row_distance;
				elseif any(vertex_sequences{i}(end) == right_side_verticies) && any(this_sequence(1) == right_side_verticies)
					expected_costs(1) = row_distance;
				else
					expected_costs(1) = vine_distance;
				end
			end
			for k=2:size_sequence
				if any(this_sequence(k-1) == left_side_verticies) && any(this_sequence(k) == left_side_verticies)
					expected_costs(k) = row_distance;
				elseif any(this_sequence(k-1) == right_side_verticies) && any(this_sequence(k) == right_side_verticies)
					expected_costs(k) = row_distance;
				else
					expected_costs(k) = vine_distance;
				end
			end
			min_costs = cost_min_frac * expected_costs;
			dist_type = repmat({cost_dist_type}, size_sequence, 1);
			if strcmp(cost_dist_type, 'Exponential')
				parameters = (1 - cost_min_frac) * ones(size_sequence, 1);
			else
				error('Distribution type not supported');
			end
			vertex_clusters{i, j, 2} = [num2cell(min_costs), dist_type, num2cell(parameters)];
			%% Create combined distribution
			if strcmp(cost_dist_type, 'Exponential')
				vertex_clusters{i, j, 3} = {sum(min_costs), 'Gamma', sum(expected_costs), 1 - cost_min_frac};
			else
				error('Distribution type not supported.');
			end
		end
	end

end


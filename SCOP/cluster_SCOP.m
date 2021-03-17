function [new_vertex_rewards, new_vertex_clusters] = cluster_SCOP(vertex_cluster_rewards, vertex_clusters, cluster_size)
%CLUSTER_SCOP Clusters verticies of a SCOP into groups
%
%	Version: 1.0
%	Date: 30/10/19
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function takes a Stochastic Cost Orienteering Problem (SCOP) and groups vertex clusters with adjacent indexes into larger clusters. Note that a cluster can be a single vertex. First used in https://ieeexplore.ieee.org/abstract/document/9340899
%	Assumptions:
%
%	Inputs:
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
%		cluster_size: the maximum number of vertex clusters to group together into a single vertex cluster
%	Outputs:
%		new_vertex_rewards: the reward associated with visiting each new vertex cluster
%		new_vertex_clusters: an upper triangular cell matrix containing information about the new vertex clusters in the same format as the input variable vertex_clusters

	%% Check inputs
	if cluster_size < 2
		new_vertex_rewards = vertex_cluster_rewards;
		new_vertex_clusters = vertex_clusters;
		return
	end
	%% Determine new vertex clusters
	n_new_vertexes = ceil((length(vertex_cluster_rewards) - 2) / cluster_size) + 2;
	if n_new_vertexes < 3
		n_new_vertexes = 3;
	end
	index = unique([1, 1+cluster_size:cluster_size:length(vertex_cluster_rewards), length(vertex_cluster_rewards)]);
	if length(index) < n_new_vertexes
		index = [index(1:end-1), index(end)-1, index(end)];
	end
	vertex_sequences = {};
	new_vertex_rewards = [];
	i = 0;
	for j=1:n_new_vertexes
		i = i + 1;
		this_sequence = [i];
		this_reward = vertex_cluster_rewards(i);
		if j~=1 & j~=n_new_vertexes
			for k=2:cluster_size
				if i+1 >= length(vertex_cluster_rewards)
					break;
				end
				i = i + 1;
				this_sequence = [this_sequence; i];
				this_reward = this_reward + vertex_cluster_rewards(i);
			end
		end
		vertex_sequences{j} = this_sequence;
		new_vertex_rewards(j) = this_reward;
	end
	%% Create new vertex clusters
	new_vertex_clusters = cell(n_new_vertexes, n_new_vertexes, 3);
	for i=1:n_new_vertexes
		for j=i+1:n_new_vertexes
			%% Build sequence of vertices
			this_vs = [vertex_clusters{vertex_sequences{i}(end), vertex_sequences{j}(1), 1}];
			this_vs_dist = [vertex_clusters{vertex_sequences{i}(end), vertex_sequences{j}(1), 2}];
			this_dist = vertex_clusters{vertex_sequences{i}(end), vertex_sequences{j}(1), 3};
			if ~strcmp(this_dist{2}, 'Gamma')
				error('Distribution type not supported.');
			end
			for k=2:length(vertex_sequences{j})
				this_vs = [this_vs; vertex_clusters{vertex_sequences{j}(k-1),vertex_sequences{j}(k), 1}];
				temp = vertex_clusters{vertex_sequences{j}(k-1),vertex_sequences{j}(k), 2};
				temp1 = [this_vs_dist(:, 1); temp(:, 1)];
				temp2 = [this_vs_dist(:, 2); temp(:, 2)];
				temp3 = [this_vs_dist(:, 3); temp(:, 3)];
				this_vs_dist = [temp1, temp2, temp3];
				%% Create combined distribution
				next_dist = vertex_clusters{vertex_sequences{j}(k-1),vertex_sequences{j}(k), 3};
				if strcmp(next_dist{2}, 'Gamma')
					this_dist{1} = this_dist{1} + next_dist{1};
					this_dist{3} = this_dist{3} + next_dist{3};
					if ~this_dist{4} == next_dist{4}
						error('Gamma distributions with different shape parameters not supported.');
					end
				else
					error('Distribution type not supported.');
				end
			end
			new_vertex_clusters{i, j, 1} = this_vs;
			new_vertex_clusters{i, j, 2} = this_vs_dist;
			new_vertex_clusters{i, j, 3} = this_dist;
		end
	end

end


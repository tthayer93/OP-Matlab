function [ tour_best, reward_best, cost_best ] = S_Algorithm_heuristic( edge_list, vertex_rewards, budget, beginning, ending, n_paths, alpha, r, allow_revisit )
%S_ALGORITHM_HEURISTIC Solves the OP using a stochastic heuristic
%
%	Version: 1.1
%	Date: 03/2/21
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves the Orienteering Problem (OP) using the Stochastic-Algorithm heuristic proposed by Tsiligirides1984. http://www.jstor.org/stable/2582629?origin=JSTOR-pdf
%	Inputs:
%		edge_list: a list containing the set of directed edges in the graph, with the following format
%			initial vertex | destination vertex | distance
%		vertex_rewards: a vector containing the reward for each vertex
%		budget: a scalar bound on the total distance of the computed route
%		beginning: the first vertex in the route
%		ending: the last vertex in the route
%		n_paths: the number of paths to build using the heuristic before returing the best, default = 3000
%		alpha: a weight factor parameter, default = 1.0
%		r: a power factor parameter, default = 4.0
%		allow_revisit: a binary variable that controls if edges can be reused. Useful if the graph is not complete. default = 0
%	Outputs:
%		tour_best: a vector containing the sequence of verticies visited by the best found path
%		reward_best: the total reward collected by the computed path
%		cost_best: the total cost associated with the computed path


    %% Initialize
	reward_list = [[1:length(vertex_rewards)]', vertex_rewards(:)]; %quick fix
	dist_list = edge_list(:, 3); %quick fix 2
    nVertex = max(reward_list(:, 1));
    if nargin < 10
        allow_revisit = 0;
        if nargin < 9
            r = 4.0;
            if nargin < 8
                alpha = 1.0;
                if nargin < 7
                    %n_paths = nVertex;
                    n_paths = 3000;
                end
            end
        end
    end
    paths = {};
    costs = [];
    total_rewards = [];
    adjacency_matrix = sparse(edge_list(:, 1), edge_list(:, 2), dist_list);
    end_dist = [];
    end_path = {};
    end_pred = [];
    for i=1:nVertex
        [end_dist(i), end_path{i}, end_pred{i}] = graphshortestpath(adjacency_matrix, i, ending);
    end
    
    %% Create Paths
    for i=1:n_paths
        current = beginning;
        tour = [beginning];
        used_budget = 0;
        rewards = reward_list(:, 2);
        total_reward = rewards(reward_list(:, 1) == beginning);
        rewards(reward_list(:, 1) == beginning) = 0;
        k=0;
        while ((used_budget + end_dist(current)) < budget) && (tour(end) ~= ending)
            lambda = 0;
            indexs = find(edge_list(:,1) == current);
            candidates = edge_list(indexs, 2);
            for j=1:length(indexs)
                lambda = lambda + rewards(reward_list(:, 1) == candidates(j))/dist_list(indexs(j));
            end
            A = zeros(nVertex, 1);
            for j=1:length(indexs)
                if ~ismember(candidates(j), tour) || allow_revisit
                    end_distance = end_dist(candidates(j));
                    E = alpha * max(budget - used_budget - dist_list(indexs(j)) - end_distance, 0) * lambda;
                    A(candidates(j)) = ((rewards(reward_list(:, 1) == candidates(j)) + E)/dist_list(indexs(j)))^r;
                end
            end
            if lambda == 0
                tour = [tour, end_path{current}(2:end)];
                total_reward = total_reward + sum(rewards(any(reward_list(:, 1) == end_path{current}(2:end), 2)));
                used_budget = used_budget + end_dist(current);
                break;
            end
            next = sum(rand >= cumsum([0; A/sum(A)]));
			index1 = edge_list(:, 1) == current;
			index2 = edge_list(:, 2) == next;
			index = index1 & index2;
			if (used_budget + dist_list(find(index, 1)) + end_dist(next)) > budget
				tour = [tour, end_path{current}(2:end)];
				total_reward = total_reward + sum(rewards(any(reward_list(:, 1) == end_path{current}(2:end), 2)));
				used_budget = used_budget + end_dist(current);
				break;
			end
			used_budget = used_budget + dist_list(find(index, 1));
			current = next;
			tour = [tour, next];
			total_reward = total_reward + rewards(reward_list(:, 1) == next);
			rewards(reward_list(:, 1) == next) = 0;
        end
        if tour(end) == ending && used_budget <= budget
            paths{i} = tour;
            costs(i) = used_budget;
            total_rewards(i) = total_reward;
        end
    end
    
    %% Find Best Solution
    tour_best = [];
    x_best = [];
    reward_best = 0;
    cost_best = 0;
    if ~isempty(total_rewards)
        [reward_best, index] = max(total_rewards);
        cost_best = costs(index);
        tour_best = paths{index};
        x_best = zeros(size(dist_list));
        for i=2:length(tour_best)
            indexs1 = edge_list(:,1) == tour_best(i-1);
            indexs2 = edge_list(:,2) == tour_best(i);
            x_best = x_best | (indexs1 & indexs2);
        end
        %tour_best = tour_best(2:end);
    end
end


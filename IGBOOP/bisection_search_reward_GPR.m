function [total_cost, total_reward, tour, avoidance_map, sampling_reward, alpha] = bisection_search_reward_GPR(vine_distance, row_distance, reward_map, beginning, ending, budget, min_rewards, sampling_map, type, avoidance_map)
%BISECTION_SEARCH_SAMPLING_GPR Solves the BOOP with reward OC on an IG
%
%	Version: 1.0
%	Date: 03/07/2019
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves the Irrigation Graph (IG) Objective Constraint (OC) Bi-Objective Orienteering Problem (BOOP) using a bisection search on the Greedy Partial Row (GPR) heuristic with minimum reward constraint, presented in https://ieeexplore.ieee.org/abstract/document/8842839
%	Assumptions:
%		The vineyard is rectangular, such that every row has the same number of vines within it.
%		There are no missing vines within the vineyard (the reward for such vines can be set to 0).
%		Vine rows are equally spaced, and each vine within a row is equally spaced.
%		The agent is allowed to turn around within a row, to access only some of the rewards within it
%		The vine_distance and row_distance are equal to 1 (if this is not the case, using avoidance_map as input or output will not be accurate)
%		Only one robot is allowed within a row at a time
%	Inputs:
%		vine_distance: The distance between each vine in the rows, which is used as the movement cost between them
%		row_distance: The distance between each row of vines, which is used as the movement cost between them
%		reward_map: A matrix of size num_rows*num_vines_per_row containing the reward of each vertex
%		beginning: The vertex at which the tour begins, which must be on the left side (column 1) of the vineyard (reward_map)
%		ending: The vertex at which the tour ends, which must be on the left side (column 1) of the vineyard (reward_map)
%		budget: The max allowable movement cost for the tour
%		min_samples: The minimum constraint for sampling rewards of the tour
%		sampling_map: A matrix of size num_rows*num_vines_per_row containing the sampling reward of each vertex
%		type: An integer defining which IGOCBOOP method to use in the bisection search
%			type = 1: Run GPR for sampling, save frac of budget for rewards, rerun GPR for rewards, binary search to find appropriate fraction of budget to save
%			type = 2; Run GPR for rewards, save frac of budget for sampling, use knapsack for sampling, binary search to find appropriate fraction of budget to save
%			type = 3; Run GPR with combined reward/sampling map using weight alpha, binary search to find appropriate alpha
%			type = 4; Run GPR for rewards, save frac of budget for sampling, add sampling points using knapsack, use up last of budget with GPR for rewards
%		avoidance_map: A cell matrix of size num_rows*num_vines_per_row, where each cell contains the times that the corresponding vertex is visited by other agents
%	Outputs:
%		total_cost: The total movement cost for the computed tour
%		total_reward: The total reward collected for the computed tour
%		tour: A sequence of vertices describing the tour computed by GPR, from beginning vertex to ending vertex
%		avoidance_map: A cell matrix of size num_rows*num_vines_per_row, where each cell contains the times that the corresponding vertex is visited by any agent
%		sampling_reward: The total sampling reward collected for the computed tour
%		alpha: A value between 0 and 1 describing the resource split between the two objectives, with reward receiving 100*alpha percent of resources and sampling receiving 100*(1-alpha) percent of resources

    rewards =  reshape(transpose(reward_map),[1, size(reward_map, 1)*size(reward_map, 2)]);
    samples = reshape(transpose(sampling_map),[1, size(sampling_map, 1)*size(sampling_map, 2)]);
    if ~exist('avoidance_map')
        avoidance_map = cell(size(reward_map));
    end
%     if type == 1
%         [total_cost, total_reward, tour, avoidance_map, sampling_reward] = dual_budget_GPR(vine_distance, row_distance, reward_map, beginning, ending, 0, budget, sampling_map, avoidance_map);
%     elseif type == 2
%         [total_cost, total_reward, tour, avoidance_map, sampling_reward] = knapsack_sampling_GPR(vine_distance, row_distance, reward_map, beginning, ending, 0, budget, sampling_map);
%     elseif type == 3
%         [total_cost, total_reward, tour, avoidance_map] = team_GPR_series(vine_distance, row_distance, sampling_map, beginning, ending, budget, avoidance_map);
%         total_reward = sum(rewards(unique(tour{1})));
%         sampling_reward = sum(samples(unique(tour{1})));
%     elseif type == 4
%         [total_cost, total_reward, tour, avoidance_map, sampling_reward] = knapsack_sampling_GPR_extended(vine_distance, row_distance, reward_map, beginning, ending, 0, budget, sampling_map);
%     end
%     if sampling_reward <= min_samples
%         alpha = 0;
%         return;
%     end
    total_cost = 0;
    total_reward = 0;
    tour = [];
    sampling_reward = 0;
    if type == 1
        [this_total_cost, this_total_reward, this_tour, this_avoidance_map, this_sampling_reward] = dual_budget_GPR(vine_distance, row_distance, reward_map, beginning, ending, 0, budget, sampling_map, avoidance_map);
    elseif type == 2
        [this_total_cost, this_total_reward, this_tour, this_avoidance_map, this_sampling_reward] = knapsack_sampling_GPR(vine_distance, row_distance, reward_map, beginning, ending, 0, budget, sampling_map);
    elseif type == 3
        [this_total_cost, this_total_reward, this_tour, this_avoidance_map] = team_GPR_series(vine_distance, row_distance, reward_map, beginning, ending, budget, avoidance_map);
        this_total_reward = sum(rewards(unique(this_tour{1})));
        this_sampling_reward = sum(samples(unique(this_tour{1})));
    elseif type == 4
        [this_total_cost, this_total_reward, this_tour, this_avoidance_map, this_sampling_reward] = knapsack_sampling_GPR_extended(vine_distance, row_distance, reward_map, beginning, ending, 0, budget, sampling_map);
    end
    if this_total_reward >= min_rewards
        total_cost = this_total_cost;
        total_reward = this_total_reward;
        tour = this_tour;
        avoidance_map = this_avoidance_map;
        sampling_reward = this_sampling_reward;
        alpha = 0;
        %return;
    end
    high_alpha = 1;
    low_alpha = 0;
    alpha = 0.5;
    while ((sampling_reward ~= this_sampling_reward) || (total_reward ~= this_total_reward)) && (abs(high_alpha - low_alpha) >= 10^-4)
        this_reward_budget = alpha*budget;
        this_sampling_budget = (1-alpha)*budget;
        if type == 1
            [this_total_cost, this_total_reward, this_tour, this_avoidance_map, this_sampling_reward] = dual_budget_GPR(vine_distance, row_distance, reward_map, beginning, ending, this_reward_budget, this_sampling_budget, sampling_map, avoidance_map);
        elseif type == 2
            [this_total_cost, this_total_reward, this_tour, this_avoidance_map, this_sampling_reward] = knapsack_sampling_GPR(vine_distance, row_distance, reward_map, beginning, ending, this_reward_budget, this_sampling_budget, sampling_map);
        elseif type == 3
            [this_total_cost, this_total_reward, this_tour, this_avoidance_map] = team_GPR_series(vine_distance, row_distance, alpha*reward_map./sum(sum(reward_map)) + (1-alpha)*sampling_map./sum(sum(sampling_map)), beginning, ending, budget, avoidance_map);
            this_total_reward = sum(rewards(unique(this_tour{1})));
            this_sampling_reward = sum(samples(unique(this_tour{1})));
        elseif type == 4
            [this_total_cost, this_total_reward, this_tour, this_avoidance_map, this_sampling_reward] = knapsack_sampling_GPR_extended(vine_distance, row_distance, reward_map, beginning, ending, this_reward_budget, this_sampling_budget, sampling_map);
        end
        if (this_total_reward >= min_rewards) && (sum(sum(sampling_map)) == this_sampling_reward)
            total_cost = this_total_cost;
            total_reward = this_total_reward;
            tour = this_tour;
            avoidance_map = this_avoidance_map;
            sampling_reward = this_sampling_reward;
            low_alpha = alpha;
            alpha = (high_alpha + low_alpha) / 2;
        elseif (this_total_reward >= min_rewards)
            total_cost = this_total_cost;
            total_reward = this_total_reward;
            tour = this_tour;
            avoidance_map = this_avoidance_map;
            sampling_reward = this_sampling_reward;
            high_alpha = alpha;
            alpha = (high_alpha + low_alpha) / 2;
        else
            high_alpha = alpha;
            alpha = (high_alpha + low_alpha) / 2;
        end
    end
end


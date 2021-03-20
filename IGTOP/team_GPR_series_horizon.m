function [ total_cost, total_reward, tours, avoidance_map ] = team_GPR_series_horizon( vine_distance, row_distance, reward_map, beginning, ending, budgets, horizon, best_of, num_moves, deeper, randomized )
%TEAM_GPR_SERIES_HORIZON Builds an orienteering tour on an IG
%
%	Version: 1.0
%	Date: 05/29/2019
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves the Irrigation Graph Orienteering Problem (IGOP) using the Greedy Partial Row (GPR) heuristic, using a greedy rolling-horizon and randomization, where decisions are full movements across the horizon.
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
%		budgets: An array containing the max allowable movement cost for the tour of each robot
%		horizon: The horizon distance at which the greedy algorithm looks
%		best_of: The number of possible moves to look at for each horizon level
%		num_moves: The number of moves for each decision, following the best found set according to the horizon depth and heuristic
%		deeper: The distance to explore past the horizion
%		randomized: The number of random movements to explore at each horizon level
%	Outputs:
%		total_cost: An array containing the total movement cost for the computed tour of each robot
%		total_reward: An array containing the total reward collected for the computed tour of each robot
%		tours: A cell array containing the computed tours for each robot, with individual cells containing a single robot's tour
%		avoidance_map: A cell matrix of size num_rows*num_vines_per_row, where each cell contains the times that the corresponding vertex is visited by any agent

    reward_map_rows = size(reward_map, 1);
    reward_map_cols = size(reward_map, 2);
    total_cost = zeros(1, length(budgets));
    rewards1 = zeros(1, length(budgets));
    tours = cell(1, length(budgets));
    reward_map1 = reward_map;
    avoidance_map = cell(size(reward_map1));
    for i=1:length(budgets)
        %[total_cost1(i), rewards1(i), tours1{i}, avoidance_map] = greedy_partial_row_avoidance(vine_distance, row_distance, reward_map1, beginning, ending, budgets(i), avoidance_map);
        [total_cost(i), rewards1(i), tours{i}, avoidance_map] = greedy_partial_row_horizon_fullmove_randomized( vine_distance, row_distance, reward_map1, beginning, ending, budgets(i), horizon, best_of, num_moves, deeper, randomized, avoidance_map);
        for j=1:length(tours{i})
            if tours{i}(j) > 0
                this_row = ceil(tours{i}(j)/reward_map_cols);
                this_col = rem(tours{i}(j) - 1, reward_map_cols) + 1;
                reward_map1(this_row, this_col) = 0;
            end
        end
    end
    total_reward = sum(rewards1);
    total_cost = sum(total_cost);

end
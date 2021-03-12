function [ total_cost, total_reward, tours ] = team_GPR_series( vine_distance, row_distance, reward_map, beginning, ending, budgets )
%TEAM_GPR_SERIES Solve the IGTOP using the GPR heuristic in series
%
%	Version: 1.0
%	Date: 03/17/2018
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves the Irrigation Graph Team Orienteering Problem (IGOP) using the Greedy Partial Row (GPR) heuristic in series, presented in https://ieeexplore.ieee.org/abstract/document/9062300
%	Assumptions:
%		The vineyard is rectangular, such that every row has the same number of vines within it.
%		There are no missing vines within the vineyard (the reward for such vines can be set to 0).
%		Vine rows are equally spaced, and each vine within a row is equally spaced.
%		The agent is allowed to turn around within a row, to access only some of the rewards within it
%		Only one robot is allowed within a row at a time
%		Each robot is to receive its own section of the vineyard to devote all of its resources to
%	Inputs:
%		vine_distance: The distance between each vine in the rows, which is used as the movement cost between them
%		row_distance: The distance between each row of vines, which is used as the movement cost between them
%		reward_map: A matrix of size num_rows*num_vines_per_row containing the reward of each vertex
%		beginning: The vertex at which the tour begins, which must be on the left side (column 1) of the vineyard (reward_map)
%		ending: The vertex at which the tour ends, which must be on the left side (column 1) of the vineyard (reward_map)
%		budgets: An array containing the max allowable movement cost for the tour of each robot
%	Outputs:
%		total_cost: An array containing the total movement cost for the computed tour of each robot
%		total_reward: An array containing the total reward collected for the computed tour of each robot
%		tours: A cell array containing the computed tours for each robot, with individual cells containing a single robot's tour

	%% Initialize
    reward_map_rows = size(reward_map, 1);
    reward_map_cols = size(reward_map, 2);
    total_cost = zeros(1, length(budgets));
    rewards1 = zeros(1, length(budgets));
    tours = cell(1, length(budgets));
    reward_map1 = reward_map;
    avoidance_map = cell(size(reward_map1));
    %% Loop for each robot
	for i=1:length(budgets)
        [total_cost(i), rewards1(i), tours{i}, avoidance_map] = greedy_partial_row_avoidance(vine_distance, row_distance, reward_map1, beginning, ending, budgets(i), avoidance_map);
        %[total_cost1(i), rewards1(i), tours1{i}, avoidance_map] = team_GPR_parallel(vine_distance, row_distance, reward_map1, beginning, ending, budgets(i), avoidance_map);
        %tours1{i} = tours1{i}{:};
        for j=1:length(tours{i})
            if tours{i}(j) > 0
                this_row = ceil(tours{i}(j)/reward_map_cols);
                this_col = rem(tours{i}(j) - 1, reward_map_cols) + 1;
                reward_map1(this_row, this_col) = 0;
            end
        end
    end
    total_reward = sum(rewards1);

end
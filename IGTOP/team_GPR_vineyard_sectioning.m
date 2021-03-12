function [ total_cost, total_reward, tours ] = team_GPR_vineyard_sectioning( vine_distance, row_distance, reward_map, beginning, ending, budgets )
%TEAM_GPR_VINEYARD_SECTIONING Solve the IGTOP using the vineyard sectioning approach
%
%	Version: 1.0
%	Date: 03/17/2018
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves the Irrigation Graph Team Orienteering Problem (IGOP) using the vineyard sectioning approach with the Greedy Partial Row (GPR) heuristic, presented in https://ieeexplore.ieee.org/abstract/document/9062300
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
    total_budget = sum(budgets);
    total_reward_map = sum(sum(reward_map));
    total_reward_rows = sum(reward_map, 2);
    total_cost = zeros(1, length(budgets));
    rewards = zeros(1, length(budgets));
    tours = cell(1, length(budgets));
	%% Loop for each robot
	current_row = 1;
    for i=1:length(budgets)
        percent_budget = budgets(i)/total_budget;
        k = current_row;
        temp = 0;
        while (temp/total_reward_map < percent_budget) && (k <= length(total_reward_rows))
            temp = temp + total_reward_rows(k);
            k = k + 1;
        end
        new_reward_map = reward_map;
        new_reward_map(1:(current_row-1), :) = 0;
        new_reward_map(k:end, :) = 0;
        %[total_cost(i), rewards(i), tours{i}] = team_GPR_parallel(vine_distance, row_distance, new_reward_map, beginning, ending, budgets(i));
        [total_cost(i), rewards(i), tours{i}] = greedy_partial_row_avoidance(vine_distance, row_distance, new_reward_map, beginning, ending, budgets(i));
        current_row = k;
    end
    total_reward = sum(rewards);
end


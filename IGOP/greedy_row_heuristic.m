function [ cost, reward, real_tour ] = greedy_row_heuristic( row_distance, vine_distance, reward_map, beginning, ending, budget )
%GREEDY_ROW_HEURISTIC Builds an orienteering tour on an IG
%
%	Version: 1.0
%	Date: 08/15/2019
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves the Irrigation Graph Orienteering Problem (IGOP) using the Greedy Row heuristic, presented in https://ieeexplore.ieee.org/abstract/document/8461242
%	Assumptions:
%		The vineyard is rectangular, such that every row has the same number of vines within it.
%		There are no missing vines within the vineyard (the reward for such vines can be set to 0).
%		Vine rows are equally spaced, and each vine within a row is equally spaced.
%		The agent must complete a row without turning around
%	Inputs:
%		row_distance: The distance between each row of vines, which is used as the movement cost between them
%		vine_distance: The distance between each vine in the rows, which is used as the movement cost between them
%		reward_map: A matrix of size num_rows*num_vines_per_row containing the reward of each vertex
%		beginning: The vertex at which the tour begins, which must be on the left side (column 1) of the vineyard (reward_map)
%		ending: The vertex at which the tour ends, which must be on the left side (column 1) of the vineyard (reward_map)
%		budget: The max allowable movement cost for the tour
%	Outputs:
%		total_cost: The total movement cost for the computed tour
%		total_reward: The total reward collected for the computed tour
%		tour: A sequence of vertices describing the tour computed by GPR, from beginning vertex to ending vertex

    %% Compile rewards and costs for each row
    beginning_row = ceil(beginning/size(reward_map,2));
    ending_row = ceil(ending/size(reward_map,2));
    num_vines_per_row = size(reward_map, 2);
    for i=1:size(reward_map,1)
        compiled_rewards(i) = sum(reward_map(i,:));
    end
    compiled_row_cost = num_vines_per_row * vine_distance;
    [rewards, index] = sort(compiled_rewards,'descend');
    rewards = [rewards', index'];
    %% Build the tour using row numbers
	current_row = beginning_row;
    current_side = 1;
    total_cost = 0;
    total_reward = 0;
    tour = [];
    while length(rewards) > 0
        if current_side == 1
            loop_cost = row_distance*abs(current_row - rewards(1, 2)) + 2*compiled_row_cost + row_distance*abs(rewards(1,2) - ending_row);
            if loop_cost <= (budget - total_cost)
                total_cost = total_cost + row_distance*abs(current_row - rewards(1, 2)) + compiled_row_cost;
                total_reward = total_reward + rewards(1, 1);
                tour = [tour; rewards(1, 2)];
                current_row = rewards(1, 2);
                current_side = 2;
                rewards(1, :) = [];
            else
                rewards(1, :) = [];
            end
        elseif current_side == 2
            loop_cost = row_distance*abs(current_row - rewards(1, 2)) + compiled_row_cost + row_distance*abs(rewards(1,2) - ending_row);
            if loop_cost <= (budget - total_cost)
                total_cost = total_cost + row_distance*abs(current_row - rewards(1, 2)) + compiled_row_cost;
                total_reward = total_reward + rewards(1, 1);
                tour = [tour; rewards(1, 2)];
                current_row = rewards(1, 2);
                current_side = 1;
                rewards(1, :) = [];
            else
                rewards(1, :) = [];
            end
        end
    end
    if current_side == 2 % Will only be at side 2 if the row needed to reach the ending has already been traversed
        total_cost = total_cost + row_distance*abs(current_row - ending_row) + compiled_row_cost;
    end
    if current_side == 1
        total_cost = total_cost + row_distance*abs(current_row - ending_row);
	end
	%% Build the complete tour using the row numbers calculated earlier
    now = beginning;
    reward = 0;
    real_tour = [];
    cost = 0;
    for i=1:length(tour)
        now_row = ceil(now/size(reward_map,2));
        now_col = rem((now-1),size(reward_map,2)) + 1;
        if now_row < tour(i)
            for j=now_row:(tour(i)-1)
                real_tour = [real_tour; (j-1)*size(reward_map,2) + now_col];
                reward = reward + reward_map(j, now_col);
                reward_map(j, now_col) = 0;
                cost = cost + row_distance;
            end
        else
            for j=now_row:-1:(tour(i)+1)
                real_tour = [real_tour; (j-1)*size(reward_map,2) + now_col];
                reward = reward + reward_map(j, now_col);
                reward_map(j, now_col) = 0;
                cost = cost + row_distance;
            end
        end
        if now_col == 1
            for j=now_col:(size(reward_map,2)-1)
                real_tour = [real_tour; (tour(i)-1)*size(reward_map,2) + j];
                reward = reward + reward_map(tour(i), j);
                reward_map(tour(i), j) = 0;
                cost = cost + vine_distance;
            end
            now = (tour(i)-1)*size(reward_map,2) + size(reward_map,2);
        else
            for j=now_col:-1:2
                real_tour = [real_tour; (tour(i)-1)*size(reward_map,2) + j];
                reward = reward + reward_map(tour(i), j);
                reward_map(tour(i), j) = 0;
                cost = cost + vine_distance;
            end
            now = (tour(i)-1)*size(reward_map,2) + 1;
        end
    end
    now_row = ceil(now/size(reward_map,2));
    now_col = rem((now-1),size(reward_map,2)) + 1;
	%% End the tour
    if now_row < ending_row
        for j=now_row:(ending_row)
            real_tour = [real_tour; (j-1)*size(reward_map,2) + now_col];
            reward = reward + reward_map(j, now_col);
            reward_map(j, now_col) = 0;
            cost = cost + row_distance;
        end
    else
        for j=now_row:-1:(ending_row)
            real_tour = [real_tour; (j-1)*size(reward_map,2) + now_col];
            reward = reward + reward_map(j, now_col);
            reward_map(j, now_col) = 0;
            cost = cost + row_distance;
        end
    end
    now = real_tour(end);
    now_row = ceil(now/size(reward_map,2));
    now_col = rem((now-1),size(reward_map,2)) + 1;
    if now_col ~= 1
        for j=now_col-1:-1:1
            real_tour = [real_tour; (now_row-1)*size(reward_map,2) + j];
            reward = reward + reward_map(now_row, j);
            reward_map(now_row, j) = 0;
            cost = cost + vine_distance;
        end
	end
    %% Eliminate duplicates
    i = 2;
    while i <= length(tour)
        if tour(i-1) == tour(i)
            tour(i-1) = [];
        else
            i = i + 1;
        end
    end
end


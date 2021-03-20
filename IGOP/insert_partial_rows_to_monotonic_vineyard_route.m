function [ total_cost, total_reward, tour ] = insert_partial_rows_to_monotonic_vineyard_route( vine_distance, row_distance, reward_map, budget, tour, total_cost, reward_goal )
%INSERT_PARTIAL_ROWS_TO_MONOTONIC Insert partial rows to monotonic vineyard route
%
%	Version: 1.0
%	Date: 02/03/2019
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function takes a monotonic vineyard route (output of make_vineyard_route_monotonic.m with omit_partials=1) and adds subtours to it to increase the total reward to the reward goal.
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
%		tour: A sequence of vertices describing the monotonic tour, from beginning vertex to ending vertex
%		total_cost: The total movement cost for the computed tour
%		reward_goal: The goal amount of reward, usually the reward of the initial tour before it was made monotonic, or inf to maximize reward
%	Outputs:
%		total_cost: The total movement cost for the computed tour
%		total_reward: The total reward collected for the computed tour
%		tour: A sequence of vertices describing the newly computed monotonic tour, from beginning vertex to ending vertex
    
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
    
    %% Figure out which rows we can expand into
    left_expandable = zeros(size_row, 1);
    right_expandable = zeros(size_row, 1);
    for i=1:length(tour)-1
        if (which_side(i) > 0) && (which_side(i+1) > 0)
            this_row = ceil(tour(i) / size_column);
            if which_side(i) == 1
                left_expandable(this_row) = 1;
            elseif which_side(i) == 2
                right_expandable(this_row) = 1;
            end
        end
    end
    
    %% Update rewards
    transverse_reward_map = reward_map';
    reward_vector = transverse_reward_map(:);
    total_reward = sum(reward_vector(unique(tour)));
    old_reward_map = reward_map;
    for i=1:length(tour)
        this_row = ceil(tour(i) / size_column);
        this_col = mod(tour(i) - 1, size_column) + 1;
        reward_map(this_row, this_col) = 0;
    end
    
    %% Calculate costs and rewards
    cumulative_reward_1 = zeros(size(reward_map));
    cumulative_reward_2 = zeros(size(reward_map));
    cumulative_cost_1 = zeros(size(reward_map));
    cumulative_cost_2 = zeros(size(reward_map));
    for i=1:size_row
        for j=1:size_column
            if left_expandable(i) == 1
                cumulative_reward_1(i, j) = sum(reward_map(i, 1:j));
            end
            if right_expandable(i) == 1
                cumulative_reward_2(i, size_column - j + 1) = sum(reward_map(i, size_column:-1:size_column-j+1));
            end
            cumulative_cost_1(i, j) = (j - 1) * vine_distance;
            cumulative_cost_2(i, size_column - j + 1) = (j - 1) * vine_distance;
        end
    end
    
    %% Perform heuristics
    while sum(sum(cumulative_reward_1)) > 0 || sum(sum(cumulative_reward_2)) > 0
        if total_reward >= reward_goal
            break;
        end
        heuristics_left = cumulative_reward_1./(2*cumulative_cost_1);
        heuristics_right = cumulative_reward_2./(2*cumulative_cost_2);
        %heuristics_left = cumulative_reward_1;
        %heuristics_right = cumulative_reward_2;
        best_heuristic_left = max(heuristics_left(:));
        best_heuristic_right = max(heuristics_right(:));
        best_heuristic = max([best_heuristic_left, best_heuristic_right]);
        if best_heuristic == best_heuristic_left
            [best_row, best_col] = find(heuristics_left == best_heuristic, 1);
            cost = 2*cumulative_cost_1(best_row, best_col);
            current_side = 1;
            direction = 1;
        else
            [best_row, best_col] = find(heuristics_right == best_heuristic, 1);
            cost = 2*cumulative_cost_2(best_row, best_col);
            current_side = size_column;
            direction = -1;
        end
        if total_cost + cost <= budget
            idx = find(tour == (best_row-1)*size_column + current_side, 1);
            subtour = [];
            cost = 0;
            for i=current_side+direction:direction:best_col
                subtour = [subtour; (best_row-1)*size_column + i];
                cost = cost + vine_distance;
            end
            for i=best_col-direction:-direction:current_side
                subtour = [subtour; (best_row-1)*size_column + i];
                cost = cost + vine_distance;
            end
            total_cost = total_cost + cost;
            total_reward = total_reward + max([direction*cumulative_reward_1(best_row, best_col),-direction*cumulative_reward_2(best_row, best_col)]);
            tour = [tour(1:idx); subtour; tour(idx+1:end)];
        end
        if current_side == 1
            cumulative_reward_1(best_row, best_col:size_column) = 0;
        else
            cumulative_reward_2(best_row, 1:best_col) = 0;
        end
    end
    
end
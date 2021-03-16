function [total_cost, total_reward, tours, avoidance_map, sampling_reward, contains_detours] = knapsack_sampling_GPR(vine_distance, row_distance, reward_map, beginning, ending, budgets, saved_budgets, sampling_map)
%KNAPSACK_SAMPLING_GPR Solves the BOOP with DM on an IG using knapsack for sampling
%
%	Version: 1.0
%	Date: 03/14/2019
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves the Irrigation Graph (IG) Dual Maximization (DM) Bi-Objective Orienteering Problem (BOOP) using the Greedy Partial Row (GPR) heuristic and adds sample locations using the knapsack problem, presented in https://ieeexplore.ieee.org/abstract/document/8842839
%	Assumptions:
%		This function only works for the single robot case. DO NOT USE FOR MULTI-ROBOT CASE
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
%		saved_budgets: An array containing the amount of budget for each robot to save for sample collection
%		sampling_map: A matrix of size num_rows*num_vines_per_row containing the sampling reward of each vertex
%	Outputs:
%		total_cost: An array containing the total movement cost for the computed tour of each robot
%		total_reward: An array containing the total reward collected for the computed tour of each robot
%		tours: A cell array containing the computed tours for each robot, with individual cells containing a single robot's tour
%		avoidance_map: A cell matrix of size num_rows*num_vines_per_row, where each cell contains the times that the corresponding vertex is visited by any agent
%		sampling_reward: The total sampling reward collected for the computed tour
%		contains_detours: The number of sample reward collection detours included
    
    %% Get initial tour
    [total_cost, total_reward, tours, avoidance_map] = team_GPR_series(vine_distance, row_distance, reward_map, beginning, ending, budgets);
    %% Check for sampling locations already visited, collect and remove
    temp_sampling_map = sampling_map';
    temp_sampling_vector = temp_sampling_map(:);
    sampling_reward = 0;
    for i=1:size(tours{1}, 1)
        if temp_sampling_vector(tours{1}(i)) > 0
            sampling_reward = sampling_reward + temp_sampling_vector(tours{1}(i));
            temp_sampling_vector(tours{1}(i)) = 0;
        end
    end
    temp_sampling_map = reshape(temp_sampling_vector, size(temp_sampling_map));
    sampling_map = temp_sampling_map';
    %% Find costs to each possible sampling location
    size_row = size(reward_map, 1);
    size_column = size(reward_map, 2);
    possibility_list = [];
    for i=1:size(sampling_map, 1)
        for j=1:size(sampling_map, 2)
            if sampling_map(i, j) > 0
%                 temp = [];
%                 for k=1:size(tours{1}, 1)
%                     this_row = ceil(tours{1}(k) / size_column);
%                     this_col = mod(tours{1}(k) - 1, size_column) + 1;
%                     if this_row == i
%                         cost = vine_distance*2*abs(this_col - j);
%                     else
%                         cost = row_distance*2*abs(this_row - i) + vine_distance*2*min((this_col-1)+(j-1), (size_column-this_col)+(size_column-j));
%                     end
%                     temp = [temp; i, j, k, sampling_map(i, j), cost];
%                 end
%                 [temp2, idx] = min(temp(:, 5));
%                 possibility_list = [possibility_list; temp(idx, :)];
                
                this_row = ceil(tours{1} / size_column);
                this_col = mod(tours{1} - 1, size_column) + 1;
                cost = zeros(length(this_row), 1);
                cost = cost + vine_distance*2*abs(this_col - j) .* (this_row == i);
                cost = cost + (row_distance*2*abs(this_row - i) + vine_distance*2*min((this_col-1)+(j-1), (size_column-this_col)+(size_column-j))) .* (this_row ~= i);
                [temp2, idx] = min(cost);
                possibility_list = [possibility_list; i, j, idx, sampling_map(i, j), cost(idx)];
            end
        end
    end
    %% Find best combination of possibilities
    modified = 1;
    while modified == 1
        modified = 0;
        %% Perform Knapsack
        if ~isempty(possibility_list)
            saved_budgets = floor(saved_budgets);
            [opt_val, opt_choice] = knapsack01_dp(possibility_list(:, 4), possibility_list(:, 5), saved_budgets);
            %% Check for overlapping sampling locations in knapsack solution
            for i=1:length(opt_choice)
                for j=(i+1):length(opt_choice)
                    % check both are same row
                    if possibility_list(i, 1) == possibility_list(j, 1) && (i~=j)
                        % check both are chosen by knapsack
                        if (opt_choice(i) == 1) && (opt_choice(j) == 1)
                            % check both add to same place in tour
                            if tours{1}(possibility_list(i, 3)) == tours{1}(possibility_list(j, 3))
                                modified = 1;
                                % check to see which is longer
                                if possibility_list(i, 5) >= possibility_list(j, 5)
                                    possibility_list(i, 4) = possibility_list(i, 4) + possibility_list(j, 4);
                                    possibility_list(j, :) = [];
                                else
                                    possibility_list(j, 4) = possibility_list(j, 4) + possibility_list(i, 4);
                                    possibility_list(i, :) = [];
                                end
                            end
                        end
                    end
                    if modified == 1
                        break;
                    end
                end
                if modified == 1
                    break;
                end
            end
        end
    end
    %% Add sampling detours to main tour
    contains_detours = 0;
    compensator = zeros(size(tours{1}));
    for x=1:size(possibility_list, 1)
        if opt_choice(x) == 1
            contains_detours = contains_detours + 1;
            i = possibility_list(x, 1);
            j = possibility_list(x, 2);
            k = possibility_list(x, 3);
            compensation = sum(compensator(1:k));
            sample_reward = possibility_list(x, 4);
            cost = possibility_list(x, 5);
            detour = [];
            this_row = ceil(tours{1}(k+compensation) / size_column);
            this_col = mod(tours{1}(k+compensation) - 1, size_column) + 1;
            if i == this_row
                if j >= this_col
                    direction_j = 1;
                elseif j < this_col
                    direction_j = -1;
                end
                for y=this_col+direction_j:direction_j:j
                    detour = [detour; (i-1)*size_column + y];
                end
                for y=j-direction_j:-direction_j:this_col
                    detour = [detour; (i-1)*size_column + y];
                end
            else
                if i > this_row
                    direction_i = 1;
                elseif i < this_row
                    direction_i = -1;
                end
                if (this_col ~= 1) && (this_col ~= size_column)
                    left_dir_cost = vine_distance*2*(this_col-1) + row_distance*2*abs(this_row-i) + vine_distance*2*(j-1);
                    right_dir_cost = vine_distance*2*(size_column-this_col) + row_distance*2*abs(this_row-i) + vine_distance*2*(size_column-j);
                    if left_dir_cost < right_dir_cost
                        direction_j = -1;
                        side = 1;
                    elseif right_dir_cost < left_dir_cost
                        direction_j = 1;
                        side = size_column;
                    end
                end
                if this_col == 1
                    direction_j = -1;
                    side = 1;
                elseif this_col == size_column
                    direction_j = 1;
                    side = size_column;
                end
                for y=this_col+direction_j:direction_j:side
                    detour = [detour; (this_row-1)*size_column + y];
                end
                for y=this_row+direction_i:direction_i:i
                    detour = [detour; (y-1)*size_column + side];
                end
                for y=side-direction_j:-direction_j:j
                    detour = [detour; (i-1)*size_column + y];
                end
                for y=j+direction_j:direction_j:side
                    detour = [detour; (i-1)*size_column + y];
                end
                for y=i-direction_i:-direction_i:this_row
                    detour = [detour; (y-1)*size_column + side];
                end
                for y=side-direction_j:-direction_j:this_col
                    detour = [detour; (this_row-1)*size_column + y];
                end
            end
            sampling_reward = sampling_reward + sample_reward;
            total_cost = total_cost + cost;
            tours{1} = [tours{1}(1:k+compensation); detour; tours{1}(k+compensation+1:end)];
            compensator(k) = compensator(k) + length(detour);
        end
    end
    transverse_reward_map = reward_map';
    reward_vector = transverse_reward_map(:);
    total_reward = sum(reward_vector(unique(tours{1})));
end
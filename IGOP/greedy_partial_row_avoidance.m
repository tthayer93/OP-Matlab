function [ total_cost, total_reward, tour, avoidance_map ] = greedy_partial_row_avoidance( vine_distance, row_distance, reward_map, beginning, ending, budget, avoidance_map )
%GREEDY_PARTIAL_ROW_AVOIDANCE Builds an orienteering tour on an IG
%
%	Version: 1.0
%	Date: 01/24/2019
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves the Irrigation Graph Orienteering Problem (IGOP) using the Greedy Partial Row (GPR) heuristic, making use of an avoidance map to avoid other agents within the IG, presented in https://ieeexplore.ieee.org/abstract/document/9062300
%	Assumptions:
%		The vineyard is rectangular, such that every row has the same number of vines within it.
%		There are no missing vines within the vineyard (the reward for such vines can be set to 0).
%		Vine rows are equally spaced, and each vine within a row is equally spaced.
%		The agent is allowed to turn around within a row, to access only some of the rewards within it
%		The vine_distance and row_distance are equal to 1 (if this is not the case, using avoidance_map as input or output will not be accurate)
%	Inputs:
%		vine_distance: The distance between each vine in the rows, which is used as the movement cost between them
%		row_distance: The distance between each row of vines, which is used as the movement cost between them
%		reward_map: A matrix of size num_rows*num_vines_per_row containing the reward of each vertex
%		beginning: The vertex at which the tour begins, which must be on the left side (column 1) of the vineyard (reward_map)
%		ending: The vertex at which the tour ends, which must be on the left side (column 1) of the vineyard (reward_map)
%		budget: The max allowable movement cost for the tour
%		avoidance_map: A cell matrix of size num_rows*num_vines_per_row, where each cell contains the times that the corresponding vertex is visited by other agents
%	Outputs:
%		total_cost: The total movement cost for the computed tour
%		total_reward: The total reward collected for the computed tour
%		tour: A sequence of vertices describing the tour computed by GPR, from beginning vertex to ending vertex
%		avoidance_map: A cell matrix of size num_rows*num_vines_per_row, where each cell contains the times that the corresponding vertex is visited by any agent

    if(nargin == 6)
        avoidance_map = cell(size(reward_map));
    end
    initial_reward_map = reward_map;
    avoidance_map_initial = avoidance_map;
    saver = 0;
    unsatisfied = 1;
    while unsatisfied
        %% Setup initial values
		waiting = saver;
        avoidance_map = avoidance_map_initial;
        total_cost = 0;
        tour = [beginning*ones(waiting+1, 1)];
        beginning_row = ceil(beginning/size(reward_map,2));
        ending_row = ceil(ending/size(reward_map,2));
        current_row = beginning_row;
        current_side = 1;
        total_reward = reward_map(current_row, 1);
        reward_map(current_row, 1) = 0;
        cumulative_reward_1 = zeros(size(reward_map));
        cumulative_reward_2 = zeros(size(reward_map));
        cumulative_cost_1 = zeros(size(reward_map));
        cumulative_cost_2 = zeros(size(reward_map));
        size_row = size(reward_map, 1); 
        size_column = size(reward_map, 2);
        for i=1:size_row
            for j=1:size_column
                cumulative_reward_1(i, j) = sum(reward_map(i, 1:j));
                cumulative_reward_2(i, size_column - j + 1) = sum(reward_map(i, size_column:-1:size_column-j+1));
                cumulative_cost_1(i, j) = (j - 1) * vine_distance;
                cumulative_cost_2(i, size_column - j + 1) = (j - 1) * vine_distance;
            end
		end
        %% Find path
        while sum(sum(cumulative_reward_1)) > 0 || sum(sum(cumulative_reward_2)) > 0
            if current_side == 1
                distance_side = row_distance * abs(current_row - (1:size_row)');
                loop_heuristic = cumulative_reward_1(:, size_column)./(cumulative_cost_1(:, size_column) + distance_side);
                muted_loop_heuristic = cumulative_reward_1./(2 * cumulative_cost_1 + distance_side);

                unsatisfied1 = 0;
                satisfied = 0;
                idk = 0;
				%% Look for available rows
                while ~satisfied
                    [best_loop_heuristic, best_loop_index] = max(loop_heuristic);
                    row_avoidance = [avoidance_map{best_loop_index, 2:size_column-1}];
                    to_loop_cost = row_distance*abs(current_row - best_loop_index);
                    through_loop_cost = cumulative_cost_1(best_loop_index, size_column);
                    to_end_cost = row_distance*abs(best_loop_index - ending_row);
                    time_enter = total_cost + to_loop_cost + waiting;
                    time_exit = time_enter + through_loop_cost + waiting;
                    if any(ismember(time_enter:time_exit, row_avoidance))
                        if (loop_heuristic(best_loop_index) > 0) && (to_loop_cost + 2*through_loop_cost + to_end_cost + 1 <= budget - total_cost)
                            idk = loop_heuristic(best_loop_index);
                        end
                        loop_heuristic(best_loop_index) = 0;
                    else
                        satisfied = 1;
                    end
                    if best_loop_heuristic == 0
                        unsatisfied1 = 1;
                        break;
                    end
                end
                loop_cost = to_loop_cost + 2*through_loop_cost + to_end_cost;
                %% Look for available partial rows
				unsatisfied2 = 0;
                satisfied = 0;
                while ~satisfied
                    best_muted_heuristic = max(muted_loop_heuristic(:));
                    [best_muted_row, best_muted_col] = find(muted_loop_heuristic == best_muted_heuristic, 1);
                    muted_avoidance = [avoidance_map{best_muted_row, 2:best_muted_col}];
                    to_muted_cost = row_distance*abs(current_row - best_muted_row);
                    through_muted_cost = 2*cumulative_cost_1(best_muted_row, best_muted_col);
                    to_end_cost = row_distance*abs(best_muted_row - ending_row);
                    time_enter = total_cost + to_muted_cost + waiting;
                    time_exit = time_enter + through_muted_cost + waiting;
                    if any(ismember(time_enter:time_exit, muted_avoidance))
                        if (muted_loop_heuristic(best_muted_row, best_muted_col) > 0) && (to_muted_cost + through_muted_cost + to_end_cost + 1 <= budget - total_cost)
                            idk = muted_loop_heuristic(best_muted_row, best_muted_col);
                        end
                        muted_loop_heuristic(best_muted_row, best_muted_col) = 0;
                    else
                        satisfied = 1;
                    end
                    if best_muted_heuristic == 0
                        unsatisfied2 = 1;
                        break;
                    end
                end
                muted_loop_cost = to_muted_cost + through_muted_cost + to_end_cost;
				%% Wait if needed
                if unsatisfied1 && unsatisfied2
                    if idk > 0
                        %total_cost = total_cost + 1;
						waiting = waiting + 1;
                        tour = [tour; tour(end)];
                    else
                        break;
                    end
				end
				%% Execute best action
                %[best_loop_heuristic, best_loop_index] = max(loop_heuristic);
                %best_muted_heuristic = max(muted_loop_heuristic(:));
                %[best_muted_row, best_muted_col] = find(muted_loop_heuristic == best_muted_heuristic, 1);
                %loop_cost = row_distance*abs(current_row - best_loop_index) + 2*cumulative_cost_1(best_loop_index, size_column) + row_distance*abs(best_loop_index - ending_row);
                %muted_loop_cost = row_distance*abs(current_row - best_muted_row) + 2*cumulative_cost_1(best_muted_row, best_muted_col) + row_distance*abs(best_muted_row - ending_row);
                if (best_loop_heuristic >= best_muted_heuristic) && (loop_cost <= (budget - total_cost))
                    %total_cost = total_cost + row_distance*abs(current_row - best_loop_index) + cumulative_cost_1(best_loop_index, size_column);
                    %total_reward = total_reward + cumulative_reward_1(best_loop_index, size_column);
                    %tour = [tour; best_loop_index, size_column];
                    if current_row <= best_loop_index
                        merp = 1;
                    else
                        merp = -1;
                    end
                    for i=current_row+merp:merp:best_loop_index
                            total_cost = total_cost + 1;
                            tour = [tour; (i-1)*size_column+1];
                            avoidance_map{i, 1} = [avoidance_map{i, 1}, total_cost + waiting];
                            total_reward = total_reward + cumulative_reward_1(i, 1);
                            cumulative_reward_1(i, :) = max(0, cumulative_reward_1(i, :) - cumulative_reward_1(i, 1));
                            cumulative_reward_2(i, 1) = cumulative_reward_2(i, 2);
                    end
                    for i=2:size_column
                        total_cost = total_cost + 1;
                        tour = [tour; (best_loop_index-1)*size_column + i]; 
                        avoidance_map{best_loop_index, i} = [avoidance_map{best_loop_index, i}, total_cost + waiting];
                    end
                    current_row = best_loop_index;
                    current_side = 2;
                    total_reward = total_reward + cumulative_reward_1(best_loop_index, size_column);
                    cumulative_reward_1(best_loop_index, :) = 0;
                    cumulative_reward_2(best_loop_index, :) = 0;
                elseif (muted_loop_cost <= (budget - total_cost))
                    %total_cost = total_cost + row_distance*abs(current_row - best_muted_row) + 2*cumulative_cost_1(best_muted_row, best_muted_col);
    %                 accumulated_reward = cumulative_reward_1(best_muted_row, best_muted_col);
    %                 total_reward = total_reward + accumulated_reward;
                    %tour = [tour; best_muted_row, best_muted_col; best_muted_row, 1];
                    if current_row <= best_muted_row
                        merp = 1;
                    else
                        merp = -1;
                    end
                    for i=current_row+merp:merp:best_muted_row
                            total_cost = total_cost + 1;
                            tour = [tour; (i-1)*size_column+1];
                            avoidance_map{i, 1} = [avoidance_map{i, 1}, total_cost + waiting];
                            total_reward = total_reward + cumulative_reward_1(i, 1);
                            cumulative_reward_1(i, :) = max(0, cumulative_reward_1(i, :) - cumulative_reward_1(i, 1));
                            cumulative_reward_2(i, 1) = cumulative_reward_2(i, 2);
                    end
                    for i=2:best_muted_col
                        total_cost = total_cost + 1;
                        tour = [tour; (best_muted_row-1)*size_column + i];
                        avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, total_cost + waiting];
                    end
                    for i=best_muted_col-1:-1:1
                        total_cost = total_cost + 1;
                        tour = [tour; (best_muted_row-1)*size_column + i];
                        avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, total_cost + waiting];
                    end
                    accumulated_reward = cumulative_reward_1(best_muted_row, best_muted_col);
                    total_reward = total_reward + accumulated_reward;
                    current_row = best_muted_row;
                    cumulative_reward_1(best_muted_row, :) = max(0, cumulative_reward_1(best_muted_row, :) - accumulated_reward);
                    %cumulative_reward_2(best_muted_row, 1:best_muted_col) = cumulative_reward_2(best_muted_row, best_muted_col + 1);
                    cumulative_reward_2(best_muted_row, 1:best_muted_col) = cumulative_reward_2(best_muted_row, min(best_muted_col + 1, size_column));
                else
                    temp_cost = row_distance*abs(current_row - (1:size_row))' + 2*cumulative_cost_1 + row_distance*abs((1:size_row) - ending_row)';
                    reachable = (temp_cost <= (budget - total_cost));
                    cumulative_reward_1(~reachable) = 0;
                    for i=1:size_row
                        if any(~reachable(i,:))
                            cumulative_reward_2(i, :) = 0;
                        end
                    end
    %                 if loop_cost > (budget - total_cost)
    %                     cumulative_reward_1(best_loop_index, :) = 0;
    %                     cumulative_reward_2(best_loop_index, :) = 0;
    %                 end
    %                 if muted_loop_cost > (budget - total_cost)
    %                     cumulative_reward_1(best_muted_row, best_muted_col:size_column) = 0;
    %                     %cumulative_reward_2(best_muted_row, :) = 0;
    %                     cumulative_reward_2(best_muted_row, :) = max(cumulative_reward_2(best_muted_row, :) - cumulative_reward_2(best_muted_row, best_muted_col), 0);
    %                 end
                end
            elseif current_side == 2
                distance_side = row_distance * abs(current_row - (1:size_row)');
                loop_heuristic = cumulative_reward_2(:, 1)./(cumulative_cost_2(:, 1) + distance_side);
                muted_loop_heuristic = cumulative_reward_2./(2 * cumulative_cost_2 + distance_side);
				%% Look for available rows
                unsatisfied1 = 0;
                satisfied = 0;
                idk = 0;
                while ~satisfied
                    [best_loop_heuristic, best_loop_index] = max(loop_heuristic);
                    row_avoidance = [avoidance_map{best_loop_index, 2:size_column-1}];
                    to_loop_cost = row_distance*abs(current_row - best_loop_index);
                    through_loop_cost = cumulative_cost_2(best_loop_index, 1);
                    to_end_cost = row_distance*abs(best_loop_index - ending_row);
                    time_enter = total_cost + to_loop_cost + waiting;
                    time_exit = time_enter + through_loop_cost + waiting;
                    if any(ismember(time_enter:time_exit, row_avoidance))
                        if (loop_heuristic(best_loop_index) > 0) && (to_loop_cost + through_loop_cost + to_end_cost + 1 <= budget - total_cost)
                            idk = loop_heuristic(best_loop_index);
                        end
                        loop_heuristic(best_loop_index) = 0;
                    else
                        satisfied = 1;
                    end
                    if best_loop_heuristic == 0
                        unsatisfied1 = 1;
                        break;
                    end
                end
                loop_cost = to_loop_cost + through_loop_cost + to_end_cost;
                %% Look for available partial rows
				unsatisfied2 = 0;
                satisfied = 0;
                while ~satisfied
                    best_muted_heuristic = max(muted_loop_heuristic(:));
                    [best_muted_row, best_muted_col] = find(muted_loop_heuristic == best_muted_heuristic, 1);
                    muted_avoidance = [avoidance_map{best_muted_row, best_muted_col:size_column-1}];
                    to_muted_cost = row_distance*abs(current_row - best_muted_row);
                    through_muted_cost = 2*cumulative_cost_2(best_muted_row, best_muted_col);
                    to_end_cost = row_distance*abs(best_muted_row - ending_row) + cumulative_cost_2(ending_row, 1);
                    time_enter = total_cost + to_muted_cost + waiting;
                    time_exit = time_enter + through_muted_cost + waiting;
                    if any(ismember(time_enter:time_exit, muted_avoidance))
                        if (muted_loop_heuristic(best_muted_row, best_muted_col) > 0) && (to_muted_cost + through_muted_cost + to_end_cost + 1 <= budget - total_cost)
                            idk = muted_loop_heuristic(best_muted_row, best_muted_col);
                        end
                        muted_loop_heuristic(best_muted_row, best_muted_col) = 0;
                    else
                        satisfied = 1;
                    end
                    if best_muted_heuristic == 0
                        unsatisfied2 = 1;
                        break;
                    end
                end
                muted_loop_cost = to_muted_cost + through_muted_cost + to_end_cost;
                %% Wait if needed
				if unsatisfied1 && unsatisfied2
                    if idk > 0
                        %total_cost = total_cost + 1;
						waiting = waiting + 1;
                        tour = [tour; tour(end)];
                    else
                        break;
                    end
				end
				%% Execute best action
                %[best_loop_heuristic, best_loop_index] = max(loop_heuristic);
                %best_muted_heuristic = max(muted_loop_heuristic(:));
                %[best_muted_row, best_muted_col] = find(muted_loop_heuristic == best_muted_heuristic, 1);
                %loop_cost = row_distance*abs(current_row - best_loop_index) + cumulative_cost_2(best_loop_index, 1) + row_distance*abs(best_loop_index - ending_row);
                %muted_loop_cost = row_distance*abs(current_row - best_muted_row) + 2*cumulative_cost_2(best_muted_row, best_muted_col) + row_distance*abs(best_muted_row - ending_row) + cumulative_cost_2(ending_row, 1);
                if (best_loop_heuristic >= best_muted_heuristic) && (loop_cost <= (budget - total_cost))
                    %total_cost = total_cost + row_distance*abs(current_row - best_loop_index) + cumulative_cost_2(best_loop_index, 1);
                    %total_reward = total_reward + cumulative_reward_2(best_loop_index, 1);
                    %tour = [tour; best_loop_index, 1];
                    if current_row <= best_loop_index
                        merp = 1;
                    else
                        merp = -1;
                    end
                    for i=current_row+merp:merp:best_loop_index
                            total_cost = total_cost + 1;
                            tour = [tour; (i)*size_column];
                            avoidance_map{i, size_column} = [avoidance_map{i, size_column}, total_cost + waiting];
                            total_reward = total_reward + cumulative_reward_2(i, size_column);
                            cumulative_reward_2(i, :) = max(0, cumulative_reward_2(i, :) - cumulative_reward_2(i, size_column));
                            cumulative_reward_1(i, size_column) = cumulative_reward_1(i, size_column-1);
                    end
                    for i=size_column-1:-1:1
                        total_cost = total_cost + 1;
                        tour = [tour; (best_loop_index-1)*size_column + i]; 
                        avoidance_map{best_loop_index, i} = [avoidance_map{best_loop_index, i}, total_cost + waiting];
                    end
                    current_row = best_loop_index;
                    current_side = 1;
                    total_reward = total_reward + cumulative_reward_2(best_loop_index, 1);
                    cumulative_reward_2(best_loop_index, :) = 0;
                    cumulative_reward_1(best_loop_index, :) = 0;
                elseif (muted_loop_cost <= (budget - total_cost))
                    %total_cost = total_cost + row_distance*abs(current_row - best_muted_row) + 2*cumulative_cost_2(best_muted_row, best_muted_col);
    %                 accumulated_reward = cumulative_reward_2(best_muted_row, best_muted_col);
    %                 total_reward = total_reward + accumulated_reward;
                    %tour = [tour; best_muted_row, best_muted_col; best_muted_row, size_column];
                    if current_row <= best_muted_row
                        merp = 1;
                    else
                        merp = -1;
                    end
                    for i=current_row+merp:merp:best_muted_row
                            total_cost = total_cost + 1;
                            tour = [tour; (i)*size_column];
                            avoidance_map{i, size_column} = [avoidance_map{i, size_column}, total_cost + waiting];
                            total_reward = total_reward + cumulative_reward_2(i, size_column);
                            cumulative_reward_2(i, :) = max(0, cumulative_reward_2(i, :) - cumulative_reward_2(i, size_column));
                            cumulative_reward_1(i, size_column) = cumulative_reward_1(i, size_column-1);
                    end
                    for i=size_column-1:-1:best_muted_col
                        total_cost = total_cost + 1;
                        tour = [tour; (best_muted_row-1)*size_column + i];
                        avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, total_cost + waiting];
                    end
                    for i=best_muted_col+1:1:size_column
                        total_cost = total_cost + 1;
                        tour = [tour; (best_muted_row-1)*size_column + i];
                        avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, total_cost + waiting];
                    end
                    accumulated_reward = cumulative_reward_2(best_muted_row, best_muted_col);
                    total_reward = total_reward + accumulated_reward;
                    current_row = best_muted_row;
                    cumulative_reward_2(best_muted_row, :) = max(0, cumulative_reward_2(best_muted_row, :) - accumulated_reward);
                    %cumulative_reward_1(best_muted_row, size_column:-1:best_muted_col) = cumulative_reward_1(best_muted_row, best_muted_col - 1);
                    cumulative_reward_1(best_muted_row, size_column:-1:best_muted_col) = cumulative_reward_1(best_muted_row, max(best_muted_col - 1, 1));
                else
                    temp_cost = row_distance*abs(current_row - (1:size_row))' + cumulative_cost_2(:, 1) + row_distance*abs((1:size_row) - ending_row)';
                    reachable = (temp_cost <= (budget - total_cost));
                    cumulative_reward_1(~reachable, :) = 0;
                    cumulative_reward_2(~reachable, :) = 0;
                    cumulative_reward_2(best_muted_row, 1:best_muted_col) = 0;
                    %cumulative_reward_1(best_muted_row, :) = 0;
                    cumulative_reward_1(best_muted_row, :) = max(cumulative_reward_1(best_muted_row, :) - cumulative_reward_1(best_muted_row, best_muted_col), 0);
    %                 if loop_cost > (budget - total_cost)
    %                     cumulative_reward_2(best_loop_index, :) = 0;
    %                     cumulative_reward_1(best_loop_index, :) = 0;
    %                 end
    %                 if muted_loop_cost > (budget - total_cost)
    %                     cumulative_reward_2(best_muted_row, 1:best_muted_col) = 0;
    %                     %cumulative_reward_1(best_muted_row, :) = 0;
    %                     cumulative_reward_1(best_muted_row, :) = max(cumulative_reward_1(best_muted_row, :) - cumulative_reward_1(best_muted_row, best_muted_col), 0);
    %                 end
                end
            end
		end
		%% Finish tour
		if current_side == 1
            %total_cost = total_cost + row_distance*abs(current_row - ending_row);
            %tour = [tour; ending_row, 1];
            if current_row > ending_row
                derp = -1;
            else
                derp = 1;
            end
            for i=current_row+derp:derp:ending_row
                total_cost = total_cost + 1;
                tour = [tour; (i-1)*size_column+1];
                avoidance_map{i, 1} = [avoidance_map{i, 1}, total_cost + waiting];
                %total_reward = total_reward + cumulative_reward_1(i,1);
            end
            %tour = [tour; ending];
            current_row = ending_row;
        end
        unsatisfied = 0;
        if current_side == 2 % Will only be at side 2 if the row needed to reach the ending has already been traversed
            %total_cost = total_cost + row_distance*abs(current_row - ending_row) + cumulative_cost_2(current_row, 1);
            %tour = [tour; ending_row, size_column; ending_row, 1];
            while tour(end) ~= ending
                if budget - total_cost < cumulative_cost_2(1,1)
                    unsatisfied = 1;
                    break;
                end
                satisfied = 0;
                home_stretch = row_distance*abs((1:size_row)' - current_row) + cumulative_cost_2(:, 1) + row_distance*abs((1:size_row)' - ending_row);
                for l=1:size_row
                    [temp, k] = min(home_stretch);
                    row_avoidance = [avoidance_map{k, 2:size_column-1}];
                    to_loop_cost = row_distance*abs(current_row - k);
                    through_loop_cost = cumulative_cost_2(k, 1);
                    to_end_cost = row_distance*abs(k - ending_row);
                    time_enter = total_cost + to_loop_cost + waiting;
                    time_exit = time_enter + through_loop_cost + waiting;
                    if ~any(ismember(time_enter:time_exit, row_avoidance))
                        satisfied = 1;
                        break;
                    else
                        home_stretch(k) = inf;
                    end
                end
                loop_cost = to_loop_cost + through_loop_cost + to_end_cost;
                if loop_cost <= (budget - total_cost) && satisfied
                    if current_row > k
                        derp = -1;
                    else
                        derp = 1;
                    end
                    for i=current_row+derp:derp:k
                        total_cost = total_cost + 1;
                        tour = [tour; (i)*size_column];
                        avoidance_map{i, size_column} = [avoidance_map{i, size_column}, total_cost + waiting];
                        %total_reward = total_reward + cumulative_reward_1(i,1);
                    end
                    for i=size_column-1:-1:1
                        total_cost = total_cost + 1;
                        tour = [tour; (k-1)*size_column + i];
                        avoidance_map{k, i} = [avoidance_map{k, i}, total_cost + waiting];
                    end
            %         for i=current_row:derp:ending_row
            %             tour = [tour; (i)*size_column];
            %             %total_reward = total_reward + cumulative_reward_2(i,num_vines_per_row);
            %         end
                    if current_row > ending_row
                        derp = -1;
                    else
                        derp = 1;
                    end
                    for i=k+derp:derp:ending_row
                        total_cost = total_cost + 1;
                        tour = [tour; (i-1)*size_column+1];
                        avoidance_map{i, 1} = [avoidance_map{i, 1}, total_cost + waiting];
                        %total_reward = total_reward + cumulative_reward_1(i,1);
                    end
                    current_row = k;
                    %tour = [tour; ending_row*size_column; ending];
                else
                    total_cost = total_cost + 1;
                    tour = [tour; tour(end)];
                end
            end
        end
        if ~unsatisfied
            break;
        end
        saver = saver + 1;
    end
    %% Eliminate duplicates
%     i = 2;
%     while i <= length(tour)
%         if tour(i-1) == tour(i)
%             tour(i-1) = [];
%         else
%             i = i + 1;
%         end
%     end
    %% Calculate reward
    transverse_reward_map = initial_reward_map';
    reward_vector = transverse_reward_map(:);
    total_reward = sum(reward_vector(unique(tour)));    
end
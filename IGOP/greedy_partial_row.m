function [ total_cost, total_reward, tour ] = greedy_partial_row( vine_distance, row_distance, reward_map, beginning, ending, budget )
%GREEDY_PARTIAL_ROW Builds an orienteering tour on an IG
%
%	Version: 1.0
%	Date: 01/24/2019
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves the Irrigation Graph Orienteering Problem (IGOP) using the Greedy Partial Row (GPR) heuristic, presented in https://ieeexplore.ieee.org/abstract/document/8461242
%	Assumptions:
%		The vineyard is rectangular, such that every row has the same number of vines within it.
%		There are no missing vines within the vineyard (the reward for such vines can be set to 0).
%		Vine rows are equally spaced, and each vine within a row is equally spaced.
%	Inputs:
%		vine_distance: The distance between each vine in the rows, which is used as the movement cost between them
%		row_distance: The distance between each row of vines, which is used as the movement cost between them
%		reward_map: A matrix of size num_rows*num_vines_per_row containing the reward of each vertex
%		beginning: The vertex at which the tour begins, which must be on the left side (column 1) of the vineyard (reward_map)
%		ending: The vertex at which the tour ends, which must be on the left side (column 1) of the vineyard (reward_map)
%		budget: The max allowable movement cost for the tour
%	Outputs:
%		total_cost: The total movement cost for the computed tour
%		total_reward: The total reward collected for the computed tour
%		tour: A sequence of vertices describing the tour computed by GPR, from beginning vertex to ending vertex

    %% Setup initial values
    beginning_row = ceil(beginning/size(reward_map,2));
    ending_row = ceil(ending/size(reward_map,2));
    current_row = beginning_row;
    current_side = 1;
    total_cost = 0;
    total_reward = reward_map(current_row, 1);
    reward_map(current_row, 1) = 0;
    tour = [beginning];
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
            [best_loop_heuristic, best_loop_index] = max(loop_heuristic);
            best_muted_heuristic = max(muted_loop_heuristic(:));
            [best_muted_row, best_muted_col] = find(muted_loop_heuristic == best_muted_heuristic, 1);
            loop_cost = row_distance*abs(current_row - best_loop_index) + 2*cumulative_cost_1(best_loop_index, size_column) + row_distance*abs(best_loop_index - ending_row);
            muted_loop_cost = row_distance*abs(current_row - best_muted_row) + 2*cumulative_cost_1(best_muted_row, best_muted_col) + row_distance*abs(best_muted_row - ending_row);
            if (best_loop_heuristic >= best_muted_heuristic) && (loop_cost <= (budget - total_cost))
                total_cost = total_cost + row_distance*abs(current_row - best_loop_index) + cumulative_cost_1(best_loop_index, size_column);
                total_reward = total_reward + cumulative_reward_1(best_loop_index, size_column);
                %tour = [tour; best_loop_index, size_column];
                if current_row <= best_loop_index
                    merp = 1;
                else
                    merp = -1;
                end
                for i=current_row:merp:best_loop_index
                    if i~=best_loop_index
                        tour = [tour; (i-1)*size_column+1];
                        total_reward = total_reward + cumulative_reward_1(i, 1);
                        cumulative_reward_1(i, :) = max(0, cumulative_reward_1(i, :) - cumulative_reward_1(i, 1));
                        cumulative_reward_2(i, 1) = cumulative_reward_2(i, 2);
                    end
                end
                for i=1:size_column
                    tour = [tour; (best_loop_index-1)*size_column + i]; 
                end
                current_row = best_loop_index;
                current_side = 2;
                cumulative_reward_1(best_loop_index, :) = 0;
                cumulative_reward_2(best_loop_index, :) = 0;
            elseif (muted_loop_cost <= (budget - total_cost))
                total_cost = total_cost + row_distance*abs(current_row - best_muted_row) + 2*cumulative_cost_1(best_muted_row, best_muted_col);
%                 accumulated_reward = cumulative_reward_1(best_muted_row, best_muted_col);
%                 total_reward = total_reward + accumulated_reward;
                %tour = [tour; best_muted_row, best_muted_col; best_muted_row, 1];
                if current_row <= best_muted_row
                    merp = 1;
                else
                    merp = -1;
                end
                for i=current_row:merp:best_muted_row
                    if i~=best_muted_row
                        tour = [tour; (i-1)*size_column+1];
                        total_reward = total_reward + cumulative_reward_1(i, 1);
                        cumulative_reward_1(i, :) = max(0, cumulative_reward_1(i, :) - cumulative_reward_1(i, 1));
                        cumulative_reward_2(i, 1) = cumulative_reward_2(i, 2);
                    end
                end
                for i=1:best_muted_col
                    tour = [tour; (best_muted_row-1)*size_column + i]; 
                end
                for i=best_muted_col:-1:1
                    tour = [tour; (best_muted_row-1)*size_column + i];
                end
                accumulated_reward = cumulative_reward_1(best_muted_row, best_muted_col);
                total_reward = total_reward + accumulated_reward;
                current_row = best_muted_row;
                cumulative_reward_1(best_muted_row, :) = max(0, cumulative_reward_1(best_muted_row, :) - accumulated_reward);
                cumulative_reward_2(best_muted_row, 1:best_muted_col) = cumulative_reward_2(best_muted_row, best_muted_col + 1);
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
            [best_loop_heuristic, best_loop_index] = max(loop_heuristic);
            best_muted_heuristic = max(muted_loop_heuristic(:));
            [best_muted_row, best_muted_col] = find(muted_loop_heuristic == best_muted_heuristic, 1);
            loop_cost = row_distance*abs(current_row - best_loop_index) + cumulative_cost_2(best_loop_index, 1) + row_distance*abs(best_loop_index - ending_row);
            muted_loop_cost = row_distance*abs(current_row - best_muted_row) + 2*cumulative_cost_2(best_muted_row, best_muted_col) + row_distance*abs(best_muted_row - ending_row) + cumulative_cost_2(ending_row, 1);
            if (best_loop_heuristic >= best_muted_heuristic) && (loop_cost <= (budget - total_cost))
                total_cost = total_cost + row_distance*abs(current_row - best_loop_index) + cumulative_cost_2(best_loop_index, 1);
                total_reward = total_reward + cumulative_reward_2(best_loop_index, 1);
                %tour = [tour; best_loop_index, 1];
                if current_row <= best_loop_index
                    merp = 1;
                else
                    merp = -1;
                end
                for i=current_row:merp:best_loop_index
                    if i~=best_loop_index
                        tour = [tour; (i)*size_column];
                        total_reward = total_reward + cumulative_reward_2(i, size_column);
                        cumulative_reward_2(i, :) = max(0, cumulative_reward_2(i, :) - cumulative_reward_2(i, size_column));
                        cumulative_reward_1(i, size_column) = cumulative_reward_1(i, size_column-1);
                    end
                end
                for i=size_column:-1:1
                    tour = [tour; (best_loop_index-1)*size_column + i]; 
                end
                current_row = best_loop_index;
                current_side = 1;
                cumulative_reward_2(best_loop_index, :) = 0;
                cumulative_reward_1(best_loop_index, :) = 0;
            elseif (muted_loop_cost <= (budget - total_cost))
                total_cost = total_cost + row_distance*abs(current_row - best_muted_row) + 2*cumulative_cost_2(best_muted_row, best_muted_col);
%                 accumulated_reward = cumulative_reward_2(best_muted_row, best_muted_col);
%                 total_reward = total_reward + accumulated_reward;
                %tour = [tour; best_muted_row, best_muted_col; best_muted_row, size_column];
                if current_row <= best_muted_row
                    merp = 1;
                else
                    merp = -1;
                end
                for i=current_row:merp:best_muted_row
                    if i~=best_muted_row
                        tour = [tour; (i)*size_column];
                        total_reward = total_reward + cumulative_reward_2(i, size_column);
                        cumulative_reward_2(i, :) = max(0, cumulative_reward_2(i, :) - cumulative_reward_2(i, size_column));
                        cumulative_reward_1(i, size_column) = cumulative_reward_1(i, size_column-1);
                    end
                end
                for i=size_column:-1:best_muted_col
                    tour = [tour; (best_muted_row-1)*size_column + i]; 
                end
                for i=best_muted_col:1:size_column
                    tour = [tour; (best_muted_row-1)*size_column + i];
                end
                accumulated_reward = cumulative_reward_2(best_muted_row, best_muted_col);
                total_reward = total_reward + accumulated_reward;
                current_row = best_muted_row;
                cumulative_reward_2(best_muted_row, :) = max(0, cumulative_reward_2(best_muted_row, :) - accumulated_reward);
                cumulative_reward_1(best_muted_row, size_column:-1:best_muted_col) = cumulative_reward_1(best_muted_row, best_muted_col - 1);
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
    if current_side == 1
        total_cost = total_cost + row_distance*abs(current_row - ending_row);
        %tour = [tour; ending_row, 1];
        if current_row > ending_row
            derp = -1;
        else
            derp = 1;
        end
        for i=current_row:derp:ending_row
            tour = [tour; (i-1)*size_column+1];
            %total_reward = total_reward + cumulative_reward_1(i,1);
        end
        %tour = [tour; ending];
    end
    if current_side == 2 % Will only be at side 2 if the row needed to reach the ending has already been traversed
        total_cost = total_cost + row_distance*abs(current_row - ending_row) + cumulative_cost_2(current_row, 1);
        %tour = [tour; ending_row, size_column; ending_row, 1];
        if current_row > ending_row
            derp = -1;
        else
            derp = 1;
        end
        for i=size_column:-1:1
            tour = [tour; (current_row-1)*size_column + i]; 
        end
%         for i=current_row:derp:ending_row
%             tour = [tour; (i)*size_column];
%             %total_reward = total_reward + cumulative_reward_2(i,num_vines_per_row);
%         end
        for i=current_row:derp:ending_row
            tour = [tour; (i-1)*size_column+1];
            %total_reward = total_reward + cumulative_reward_1(i,1);
        end
        %tour = [tour; ending_row*size_column; ending];
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
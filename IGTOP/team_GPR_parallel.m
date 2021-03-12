function [ total_cost, total_reward, tours, avoidance_map ] = team_GPR_parallel( vine_distance, row_distance, reward_map, beginning, ending, budgets, avoidance_map_initial )
%TEAM_GPR_PARALLEL Solve the IGTOP using the Parallel GPR heuristic
%
%	Version: 1.0
%	Date: 03/17/2018
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves the Irrigation Graph Team Orienteering Problem (IGOP) using the Parallel Greedy Partial Row (GPR) heuristic, presented in https://ieeexplore.ieee.org/abstract/document/9062300
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
%		avoidance_map_initial: A cell matrix of size num_rows*num_vines_per_row, where each cell contains the times that the corresponding vertex is visited by other agents
%	Outputs:
%		total_cost: An array containing the total movement cost for the computed tour of each robot
%		total_reward: An array containing the total reward collected for the computed tour of each robot
%		tours: A cell array containing the computed tours for each robot, with individual cells containing a single robot's tour
%		avoidance_map: A cell matrix of size num_rows*num_vines_per_row, where each cell contains the times that the corresponding vertex is visited by any agent

	%% Initialize
    if(nargin == 6)
        avoidance_map_initial = cell(size(reward_map));
	end
    waiting = zeros(1, length(budgets));
    saver = zeros(1, length(budgets));
    unsatisfied = 1;
	%% Loop incase set of tours is not viable
    while unsatisfied
		%% Initialize
        avoidance_map = avoidance_map_initial;
        nAgents = length(budgets);
        beginning_row = ceil(beginning/size(reward_map, 2));
        ending_row = ceil(ending/size(reward_map, 2));
        current_row = beginning_row * ones(1, nAgents);
        current_side = ones(1, nAgents);
        total_cost = zeros(1, nAgents);
        total_reward = reward_map(current_row(1), 1);
        reward_map(current_row, 1) = 0;
        tours = cell(1, nAgents);
        for i=1:nAgents
            tours{i} = [beginning*ones(waiting(i) + 1, 1)];
        end
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
        %% Iteratively find paths
        agent_moved = ones(1, nAgents);
        action_list = [];
        black_list = [0, 0, 0, 0, 0];
        while sum(sum(cumulative_reward_1)) > 0 || sum(sum(cumulative_reward_2)) > 0
			for i=1:nAgents
               %if agent_moved(i) == 1
					%% Compute heuristic values
                    if ~isempty(action_list)
                        action_list(action_list(:, 1) == i, :) = [];
                    end
                    distance_side = row_distance * abs(current_row(i) - (1:size_row)');
                    if current_side(i) == 1
                        loop_heuristic = cumulative_reward_1(:, size_column)./(cumulative_cost_1(:, size_column) + distance_side);
                        muted_loop_heuristic = cumulative_reward_1./(2 * cumulative_cost_1 + distance_side);
                        agent_loop_cost = row_distance*abs(current_row(i) - (1:size_row))' + 2*cumulative_cost_1(:, size_column) + row_distance*abs((1:size_row) - ending_row)';
                        agent_muted_loop_cost = row_distance*abs(current_row(i) - (1:size_row))' + 2*cumulative_cost_1 + row_distance*abs((1:size_row) - ending_row)';
                    elseif current_side(i) == 2
                        loop_heuristic = cumulative_reward_2(:, 1)./(cumulative_cost_2(:, 1) + distance_side);
                        muted_loop_heuristic = cumulative_reward_2./(2 * cumulative_cost_2 + distance_side);
                        agent_loop_cost = row_distance*abs(current_row(i) - (1:size_row))' + cumulative_cost_2(:, 1) + row_distance*abs((1:size_row) - ending_row)';
                        agent_muted_loop_cost = row_distance*abs(current_row(i) - (1:size_row))' + 2*cumulative_cost_2 + row_distance*abs((1:size_row) - ending_row)' + cumulative_cost_2(ending_row, 1);
					end
					%% Find the best viable full loop movement for this agent
                    unsatisfied1 = 0;
                    satisfied = 0;
                    idk = 0;
                    while ~satisfied
                        [best_loop_heuristic, best_loop_index] = max(loop_heuristic);
                        if best_loop_heuristic == 0
                            unsatisfied1 = 1;
                            break;
                        end
                        row_avoidance = [avoidance_map{best_loop_index, 2:size_column-1}];
                        to_loop_cost = row_distance*abs(current_row(i) - best_loop_index);
                        if current_side(i) == 1
                            through_loop_cost = cumulative_cost_1(best_loop_index, size_column);
                        else
                            through_loop_cost = cumulative_cost_2(best_loop_index, 1);
                        end
                        to_end_cost = row_distance*abs(best_loop_index - ending_row);
                        time_enter = total_cost(i) + to_loop_cost + waiting(i);
                        time_exit = time_enter + through_loop_cost;
                        checked = all((black_list(:, 1:3) == [i, best_loop_index, 0]), 2);
                        if any(ismember(time_enter:time_exit, row_avoidance)) || any(checked)
                            if (loop_heuristic(best_loop_index) > 0) && (agent_loop_cost(best_loop_index) + 1 <= budgets(i) - total_cost(i) - saver(i))
                                idk = loop_heuristic(best_loop_index);
                            end
                            loop_heuristic(best_loop_index) = 0;
                        else
                            action_list = [action_list; i, best_loop_index, 0, best_loop_heuristic, agent_loop_cost(best_loop_index)];
                            satisfied = 1;
                        end
					end
					%% Find the best viable partial loop movement for this agent
                    unsatisfied2 = 0;
                    satisfied = 0;
                    while ~satisfied
                        best_muted_heuristic = max(muted_loop_heuristic(:));
                        if best_muted_heuristic == 0
                            unsatisfied2 = 1;
                            break;
                        end
                        [best_muted_row, best_muted_col] = find(muted_loop_heuristic == best_muted_heuristic, 1);
                        to_muted_cost = row_distance*abs(current_row(i) - best_muted_row);
                        if current_side(i)
                            muted_avoidance = [avoidance_map{best_muted_row, 2:best_muted_col}];
                            through_muted_cost = 2*cumulative_cost_1(best_muted_row, best_muted_col);
                            to_end_cost = row_distance*abs(best_muted_row - ending_row);
                        else
                            muted_avoidance = [avoidance_map{best_muted_row, best_muted_col:size_column-1}];
                            through_muted_cost = 2*cumulative_cost_2(best_muted_row, best_muted_col);
                            to_end_cost = row_distance*abs(best_muted_row - ending_row) + cumulative_cost_2(ending_row, 1);
                        end
                        time_enter = total_cost(i) + to_muted_cost + waiting(i);
                        time_exit = time_enter + through_muted_cost;
                        checked = all((black_list(:, 1:3) == [i, best_muted_row, best_muted_col]), 2);
                        if any(ismember(time_enter:time_exit, muted_avoidance)) || any(checked)
                            if (muted_loop_heuristic(best_muted_row, best_muted_col) > 0) && (agent_muted_loop_cost(best_muted_row, best_muted_col) + 1 <= budgets(i) - total_cost(i) - saver(i))
                                idk = muted_loop_heuristic(best_muted_row, best_muted_col);
                            end
                            muted_loop_heuristic(best_muted_row, best_muted_col) = 0;
                        else
                            action_list = [action_list; i, best_muted_row, best_muted_col, best_muted_heuristic, agent_muted_loop_cost(best_muted_row, best_muted_col)];
                            satisfied = 1;
                        end
					end
					%% Check if viable movement found
                    agent_moved(i) = 0;
                    if unsatisfied1 && unsatisfied2
                        if idk > 0
                            total_cost(i) = total_cost(i) + 1;
                            tours{i} = [tours{i}; tours{i}(end)];
                            %agent_moved(i) = 1;
                        else
                            %break;
                        end
                    end
               %end
            end
            if isempty(action_list)
                break;
            end
            %% Make best possible movement
            while sum(agent_moved) == 0 && ~isempty(action_list)
                [temp, index] = max(action_list(:, 4));
                i = action_list(index, 1);
                j = action_list(index, 2);
                k = action_list(index, 3);
                heuristic = action_list(index, 4);
                infered_cost = action_list(index, 5);
                if infered_cost <= (budgets(i) - total_cost(i) - saver(i))
                    if current_row(i) <= j
                        merp = 1;
                    else
                        merp = -1;
					end
                    if current_side(i) == 1
						%% Traverse the left side (column 1)
                        for x=current_row(i)+merp:merp:j
                            total_cost(i) = total_cost(i) + 1;
                            tours{i} = [tours{i}; (x-1)*size_column+1];
                            avoidance_map{x, 1} = [avoidance_map{x, 1}, total_cost(i) + waiting(i)];
                            total_reward = total_reward + cumulative_reward_1(x, 1);
                            cumulative_reward_1(x, :) = max(0, cumulative_reward_1(x, :) - cumulative_reward_1(x, 1));
                            cumulative_reward_2(x, 1) = cumulative_reward_2(x, 2);
                        end
                        %total_cost(i) = total_cost(i) + row_distance*abs(current_row(i) - j);
                        %if row, else muted
                        if k == 0
                            %total_cost(i) = total_cost(i) + cumulative_cost_1(j, size_column);
                            total_reward = total_reward + cumulative_reward_1(j, size_column);
                            %add new verticies to path
                            for x=2:size_column
                                total_cost(i) = total_cost(i) + 1;
                                tours{i} = [tours{i}; (j-1)*size_column + x];
                                avoidance_map{j, x} = [avoidance_map{j, x}, total_cost(i) + waiting(i)];
                            end
                            current_row(i) = j;
                            current_side(i) = 2;
                            cumulative_reward_1(j, :) = 0;
                            cumulative_reward_2(j, :) = 0;
                        else
                            %total_cost(i) = total_cost(i) + 2*cumulative_cost_1(j, k);
                            total_reward = total_reward + cumulative_reward_1(j, k);
                            %add new verticies to path
                            for x=2:k
                                total_cost(i) = total_cost(i) + 1;
                                tours{i} = [tours{i}; (j-1)*size_column + x];
                                avoidance_map{j, x} = [avoidance_map{j, x}, total_cost(i) + waiting(i)];
                            end
                            for x=k-1:-1:1
                                total_cost(i) = total_cost(i) + 1;
                                tours{i} = [tours{i}; (j-1)*size_column + x];
                                avoidance_map{j, x} = [avoidance_map{j, x}, total_cost(i) + waiting(i)];
                            end
                            current_row(i) = j;
                            cumulative_reward_1(j, :) = max(0, cumulative_reward_1(j, :) - cumulative_reward_1(j, k));
                            cumulative_reward_2(j, 1:k) = cumulative_reward_2(j, max(k+1, size_column));
						end
					else
						%% Traverse the right side (last column)
                        for x=current_row(i)+merp:merp:j
                            total_cost(i) = total_cost(i) + 1;
                            tours{i} = [tours{i}; (x)*size_column];
                            avoidance_map{x, size_column} = [avoidance_map{x, size_column}, total_cost(i) + waiting(i)];
                            total_reward = total_reward + cumulative_reward_2(x, size_column);
                            cumulative_reward_2(x, :) = max(0, cumulative_reward_2(x, :) - cumulative_reward_2(x, size_column));
                            cumulative_reward_1(x, size_column) = cumulative_reward_1(x, size_column-1);
                        end
                        %total_cost(i) = total_cost(i) + row_distance*abs(current_row(i) - j);
                        %if row, else muted
                        if k == 0
                            %total_cost(i) = total_cost(i) + cumulative_cost_2(j, 1);
                            total_reward = total_reward + cumulative_reward_2(j, 1);
                            %add new verticies to path
                            for x=size_column-1:-1:1
                                total_cost(i) = total_cost(i) + 1;
                                tours{i} = [tours{i}; (j-1)*size_column + x];
                                avoidance_map{j, x} = [avoidance_map{j, x}, total_cost(i) + waiting(i)];
                            end
                            current_row(i) = j;
                            current_side(i) = 1;
                            cumulative_reward_2(j, :) = 0;
                            cumulative_reward_1(j, :) = 0;
                        else
                            %total_cost(i) = total_cost(i) + 2*cumulative_cost_2(j, k);
                            total_reward = total_reward + cumulative_reward_2(j, k);
                            %add new verticies to path
                            for x=size_column-1:-1:k
                                total_cost(i) = total_cost(i) + 1;
                                tours{i} = [tours{i}; (j-1)*size_column + x];
                                avoidance_map{j, x} = [avoidance_map{j, x}, total_cost(i) + waiting(i)];
                            end
                            for x=k+1:size_column
                                total_cost(i) = total_cost(i) + 1;
                                tours{i} = [tours{i}; (j-1)*size_column + x];
                                avoidance_map{j, x} = [avoidance_map{j, x}, total_cost(i) + waiting(i)];
                            end
                            current_row(i) = j;
                            cumulative_reward_2(j, :) = max(0, cumulative_reward_2(j, :) - cumulative_reward_2(j, k));
                            cumulative_reward_1(j, size_column:-1:k) = cumulative_reward_1(j, max(k-1, 1));
                        end
                    end
                    agent_moved(i) = 1;
                %if not enough budget, check if reachable by any other agent
                else
                    black_list = [black_list; action_list(index, :)];
                end
                action_list(index, :) = [];
            end
            %% Check if items on blacklist are unreachable
            %action_list = [action_list; i, best_muted_row, best_muted_col, best_muted_heuristic, agent_muted_loop_cost(best_muted_row, best_muted_col)];
            reachable = zeros(size_row, size_column);
            for x=1:nAgents
                if current_side(x) == 1
                    temp_cost = row_distance*abs(current_row(x) - (1:size_row))' + 2*cumulative_cost_1 + row_distance*abs((1:size_row) - ending_row)';
                elseif current_side(x) == 2
                    temp_cost = row_distance*abs(current_row(x) - (1:size_row))' + cumulative_cost_2(:, 1) + row_distance*abs((1:size_row) - ending_row)';
                end
                this_reachable = (temp_cost <= (budgets(x) - total_cost(x)));
                reduced_blacklist_idxs = find(black_list(:,1) == x);
                reduced_blacklist = black_list(reduced_blacklist_idxs,:);
                for y=1:size(reduced_blacklist, 1)
                    if reduced_blacklist(y, 3) == 0
                        this_reachable(reduced_blacklist(y, 2), :) = 0;
                    else
                        if current_side(x) == 1
                            this_reachable(reduced_blacklist(y, 2), 1:reduced_blacklist(y, 3)) = 0;
                        elseif current_side(x) == 2
                            this_reachable(reduced_blacklist(y, 2), size_column:-1:reduced_blacklist(y, 3)) = 0;
                        end
                    end
                end
                %black_list(reduced_blacklist_idxs, :) = [];
                reachable = ( this_reachable | reachable);
            end
            cumulative_reward_1(~reachable) = 0;
            for x=1:size_row
                if any(~reachable(x, :))
                    cumulative_reward_2(x, :) = 0;
                end
            end
		end
        %% Finish paths
        unsatisfied = 0;
        for i=1:nAgents
            while tours{i}(end) < 0
                total_cost(i) = total_cost(i) - 1;
                tours{i}(end) = [];
            end
            if current_side(i) == 2
                %total_cost(i) = total_cost(i) + row_distance*abs(current_row(i) - ending_row) + cumulative_cost_2(current_row(i), 1);
                while current_side(i) ~= 1
                    if budgets(i) - total_cost(i) < cumulative_cost_2(current_row(i), 1)
                        unsatisfied = 1;
                        %waiting(i) = waiting(i) + 1;
                        saver(i) = saver(i) + 1;
                        break;
                    end
                    home_stretch = row_distance*abs((1:size_row)' - current_row(i)) + cumulative_cost_2(current_row(i), 1) + row_distance*abs((1:size_row)' - ending_row);
                    for l=1:size_row
                        [temp, k] = min(home_stretch);
                        row_avoidance = [avoidance_map{k, 2:size_column-1}];
                        to_loop_cost = row_distance*abs(current_row(i) - k);
                        through_loop_cost = cumulative_cost_2(k, 1);
                        to_end_cost = row_distance*abs(k - ending_row);
                        time_enter = total_cost(i) + to_loop_cost + waiting(i);
                        time_exit = time_enter + through_loop_cost;
                        if ~any(ismember(time_enter:time_exit, row_avoidance))
                            satisfied = 1;
                            break;
                        else
                            home_stretch(k) = inf;
                        end
                    end
                    loop_cost = to_loop_cost + through_loop_cost + to_end_cost;
                    if loop_cost <= (budgets(i) - total_cost(i)) && satisfied
                        if current_row(i) > k
                            derp = -1;
                        else
                            derp = 1;
                        end
                        for j=current_row(i)+derp:derp:k
                            total_cost(i) = total_cost(i) + 1;
                            tours{i} = [tours{i}; (j)*size_column];
                            avoidance_map{j, size_column} = [avoidance_map{j, size_column}, total_cost(i) + waiting(i)];
                        end
                        for j=size_column-1:-1:1
                            total_cost(i) = total_cost(i) + 1;
                            tours{i} = [tours{i}; (k-1)*size_column + j];
                            avoidance_map{k, j} = [avoidance_map{k, j}, total_cost(i) + waiting(i)];
                        end
                        current_side(i) = 1;
                        current_row(i) = k;
                    else
                        total_cost(i) = total_cost(i) + 1;
                        tours{i} = [tours{i}; tours{i}(end)];
                    end
                end
                if unsatisfied
                    break;
                end
    %             for j=size_column-1:-1:1
    %                 total_cost(i) = total_cost(i) + 1;
    %                 tours{i} = [tours{i}; (current_row(i)-1)*size_column + j];
    %                 avoidance_map{current_row(i), j} = [avoidance_map{current_row(i), j}, total_cost(i)];
    %             end
            end
            if current_row(i) > ending_row
                derp = -1;
            else
                derp = 1;
            end
            for j=current_row(i)+derp:derp:ending_row
                total_cost(i) = total_cost(i) + 1;
                tours{i} = [tours{i}; (j-1)*size_column+1];
                avoidance_map{j, 1} = [avoidance_map{j, 1}, total_cost(i) + waiting(i)];
            end
            %eliminate duplicates
    %         j = 2;
    %         while j <= length(tours{i})
    %             if tours{i}(j-1) == tours{i}(j)
    %                 tours{i}(j-1) = [];
    %             else
    %                 j = j + 1;
    %             end
    %         end
        end
        if ~unsatisfied
            break;
        end
    end

end


function [total_cost, total_reward, tour, avoidance_map] = greedy_partial_row_horizon_fullmove_randomized( vine_distance, row_distance, reward_map, beginning, ending, budget, horizon, best_of, num_moves, deeper, randomized, avoidance_map)
%GREEDY_PARTIAL_ROW_HORIZON_FULLMOVE_RANDOMIZED Builds an orienteering tour on an IG
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
%		budget: The max allowable movement cost for the tour
%		horizon: The horizon distance at which the greedy algorithm looks
%		best_of: The number of possible moves to look at for each horizon level
%		num_moves: The number of moves for each decision, following the best found set according to the horizon depth and heuristic
%		deeper: The distance to explore past the horizion
%		randomized: The number of random movements to explore at each horizon level
%		avoidance_map: A cell matrix of size num_rows*num_vines_per_row, where each cell contains the times that the corresponding vertex is visited by other agents
%	Outputs:
%		total_cost: The total movement cost for the computed tour
%		total_reward: The total reward collected for the computed tour
%		tour: A sequence of vertices describing the tour computed by GPR, from beginning vertex to ending vertex
%		avoidance_map: A cell matrix of size num_rows*num_vines_per_row, where each cell contains the times that the corresponding vertex is visited by any agent

    %% Initialization
    if deeper < 0
        deeper = 0;
    end
    if horizon < 1
        horizon = 1;
    end
    if num_moves > horizon
        num_moves = horizon;
    end
    if num_moves < 1
        num_moves = 1;
    end
    if best_of < 2
        best_of = 2;
    end
    old_reward_map = reward_map;
    avoidance_map_initial = avoidance_map;
    saver = 0;
    unsatisfied = 1;
    while unsatisfied
        %% Setup initial values
        avoidance_map = avoidance_map_initial;
        total_cost = 0;
        tour = [beginning*ones(saver+1, 1)];
        waiting = 0;
        beginning_row = ceil(beginning/size(reward_map, 2));
        ending_row = ceil(ending/size(reward_map, 2));
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
        original_best_of = best_of;
        while sum(sum(cumulative_reward_1)) > 0 || sum(sum(cumulative_reward_2)) > 0
            % look across horizon
            % [1 row, 2 column, 3 ending side, 4 cost, 5 reward, 6 heuristic, 7 total_cost, 8 total_reward, 9 total_heuristic, 10 prev decision]
            decision_tree = cell(horizon+deeper, 1);
            best_of = original_best_of;
            for i=1:(horizon+deeper)
                if i > horizon
                    best_of = 1;
                end
                if i == 1
                    % decision tree root
                    branch = [current_row, -1, current_side, 0, 0, 0, 0, 0, 0, 1];
                    prev_cumulative_reward_1{1} = cumulative_reward_1;
                    prev_cumulative_reward_2{1} = cumulative_reward_2;
                else
                    % pull data for decision tree
                    branch = decision_tree{i-1};
                    prev_cumulative_reward_1 = next_cumulative_reward_1;
                    prev_cumulative_reward_2 = next_cumulative_reward_2;
                    next_cumulative_reward_1 = cell(size(branch, 1), 1);
                    next_cumulative_reward_2 = cell(size(branch, 1), 1);
                end
                % iterate over decision tree branches
                % does not account for time avoidance
                % need to account for path to visited spots on future branches
                horizon_branches = [];
                for k=1:size(branch, 1)
                    this_row = branch(k, 1);
                    this_col = branch(k, 2);
                    this_side = branch(k, 3);
                    this_cost = branch(k, 4);
                    this_reward = branch(k, 5);
                    this_heuristic = branch(k, 6);
                    this_total_cost = branch(k, 7);
                    this_total_reward = branch(k, 8);
                    this_total_heuristic = branch(k, 9);
                    temp_cumulative_reward_1 = prev_cumulative_reward_1{branch(k, 10)};
                    temp_cumulative_reward_2 = prev_cumulative_reward_2{branch(k, 10)};
                    if this_col == 0
                        temp_cumulative_reward_1(this_row, :) = 0;
                        temp_cumulative_reward_2(this_row, :) = 0;
                    elseif this_col > 0
                        if this_side == 1
                            temp_cumulative_reward_1(this_row, :) = max(0, temp_cumulative_reward_1(this_row, :) - temp_cumulative_reward_1(this_row, this_col));
                            temp_cumulative_reward_2(this_row, 1:this_col) = temp_cumulative_reward_2(this_row, min(this_col+1, size_column));
                            if this_col == size_column
                                temp_cumulative_reward_2(this_row, :) = 0;
                            end
                        elseif this_side == 2
                            temp_cumulative_reward_2(this_row, :) = max(0, temp_cumulative_reward_2(this_row, :) - temp_cumulative_reward_2(this_row, this_col));
                            temp_cumulative_reward_1(this_row, size_column:-1:this_col) = temp_cumulative_reward_1(this_row, max(this_col-1, 1));
                            if this_col == 1
                                temp_cumulative_reward_1(this_row, :) = 0;
                            end
                        end
                    end
                    if this_side == 1
                        distance_side = row_distance * abs(this_row - (1:size_row)');
                        % full row
                        temp_loop_cost = (cumulative_cost_1(:, size_column) + distance_side);
                        temp_loop_reward = temp_cumulative_reward_1(:, size_column);
                        temp_loop_heuristic = temp_loop_reward ./ temp_loop_cost;
                        temp_loop_total_cost = temp_loop_cost + this_total_cost;
                        temp_loop_total_reward = temp_loop_reward + this_total_reward;
                        temp_loop_total_heuristic = temp_loop_total_reward ./ temp_loop_total_cost;
                        loop_to_end_cost = cumulative_cost_2(:, 1) + row_distance * abs((1:size_row)' - ending_row);
                        loop_feasibility = ((temp_loop_cost + loop_to_end_cost) <= (budget - this_total_cost - total_cost));
                        temp_row_branches = [(1:size_row)', zeros(size_row, 1), 2*ones(size_row, 1), temp_loop_cost, temp_loop_reward, temp_loop_heuristic, temp_loop_total_cost, temp_loop_total_reward, temp_loop_total_heuristic, k*ones(size_row, 1)];
                        temp_row_branches = temp_row_branches(loop_feasibility == 1, :);
                        % muted row
                        temp_muted_cost = (2 * cumulative_cost_1 + distance_side);
                        temp_muted_reward = temp_cumulative_reward_1;
                        muted_to_end_cost = row_distance * abs((1:size_row)' - ending_row);
                        muted_feasibility = ((temp_muted_cost + muted_to_end_cost) <= (budget - this_total_cost - total_cost));
                        rows = reshape((1:size_row)' * ones(1, size_column), [], 1);
                        cols = reshape(ones(size_row, 1) * (1:size_column), [], 1);
                        sids = reshape(1 * ones(size_row, size_column), [], 1);
                        costs = reshape(temp_muted_cost, [], 1);
                        rewards = reshape(temp_muted_reward, [], 1);
                        heuristics = rewards ./ costs;
                        t_costs = costs + this_total_cost;
                        t_rewards = rewards + this_total_reward;
                        t_heuristics = t_rewards ./ t_costs;
                        ks = reshape(k * ones(size_row, size_column), [], 1);
                        muted_feasibility = reshape(muted_feasibility, [], 1);
                        temp_muted_branches = [rows, cols, sids, costs, rewards, heuristics, t_costs, t_rewards, t_heuristics, ks];
                        temp_muted_branches = temp_muted_branches(muted_feasibility == 1, :);
                    elseif this_side == 2
                        distance_side = row_distance * abs(this_row - (1:size_row)');
                        % full row
                        temp_loop_cost = (cumulative_cost_2(:, 1) + distance_side);
                        temp_loop_reward = temp_cumulative_reward_2(:, 1);
                        temp_loop_heuristic = temp_loop_reward ./ temp_loop_cost;
                        temp_loop_total_cost = temp_loop_cost + this_total_cost;
                        temp_loop_total_reward = temp_loop_reward + this_total_reward;
                        temp_loop_total_heuristic = temp_loop_total_reward ./ temp_loop_total_cost;
                        loop_to_end_cost = row_distance * abs((1:size_row)' - ending_row);
                        loop_feasibility = ((temp_loop_cost + loop_to_end_cost) <= (budget - this_total_cost - total_cost));
                        temp_row_branches = [(1:size_row)', zeros(size_row, 1), 1*ones(size_row, 1), temp_loop_cost, temp_loop_reward, temp_loop_heuristic, temp_loop_total_cost, temp_loop_total_reward, temp_loop_total_heuristic, k*ones(size_row, 1)];
                        temp_row_branches = temp_row_branches(loop_feasibility == 1, :);
                        % muted row
                        temp_muted_cost = (2 * cumulative_cost_2 + distance_side);
                        temp_muted_reward = temp_cumulative_reward_2;
                        muted_to_end_cost = cumulative_cost_2(:, 1) + row_distance * abs((1:size_row)' - ending_row);
                        muted_feasibility = ((temp_muted_cost + muted_to_end_cost) <= (budget - this_total_cost - total_cost));
                        rows = reshape((1:size_row)' * ones(1, size_column), [], 1);
                        cols = reshape(ones(size_row, 1) * (1:size_column), [], 1);
                        sids = reshape(2 * ones(size_row, size_column), [], 1);
                        costs = reshape(temp_muted_cost, [], 1);
                        rewards = reshape(temp_muted_reward, [], 1);
                        heuristics = rewards ./ costs;
                        t_costs = costs + this_total_cost;
                        t_rewards = rewards + this_total_reward;
                        t_heuristics = t_rewards ./ t_costs;
                        ks = reshape(k * ones(size_row, size_column), [], 1);
                        muted_feasibility = reshape(muted_feasibility, [], 1);
                        temp_muted_branches = [rows, cols, sids, costs, rewards, heuristics, t_costs, t_rewards, t_heuristics, ks];
                        temp_muted_branches = temp_muted_branches(muted_feasibility == 1, :);
                    end
                    %Eliminate certain branches
                    if ~isempty(temp_row_branches)
                        temp_row_branches = temp_row_branches(temp_row_branches(:, 4) > 0, :);
                        %temp_row_branches = temp_row_branches(temp_row_branches(:, 5) > 0, :);
                    end
                    if ~isempty(temp_muted_branches)
                        temp_muted_branches = temp_muted_branches(temp_muted_branches(:, 4) > 0, :);
                        temp_muted_branches = temp_muted_branches(temp_muted_branches(:, 5) > 0, :);
                    end
                    %Find maxes
                    top_row = [];
                    if size(temp_row_branches, 1) >= 2
                        [val, idx] = max(temp_row_branches(:, 9));
                        top_row = temp_row_branches(idx, :);
                        temp_row_branches(idx, :) = [];
                    end
                    top_muted = [];
                    if size(temp_muted_branches, 1) >= 2
                        [val, idx] = max(temp_muted_branches(:, 9));
                        top_muted = temp_muted_branches(idx, :);
                        temp_muted_branches(idx, :) = [];
                    end
                    %Sort and get best_of
                    new_branches = [temp_row_branches; temp_muted_branches];
                    if randomized <= 1
                        num_to_keep = best_of - size(top_row, 1) - size(top_muted, 1);
                    else
                        num_to_keep = best_of * randomized - size(top_row, 1) - size(top_muted, 1);
                    end
                    if ~isempty(new_branches)
                        [vals, idxes] = maxk(new_branches(:, 9), num_to_keep);
                        new_branches = new_branches(idxes, :);
                    end
                    new_branches = [top_row; top_muted; new_branches];
                    %check for avoidance
                    keeper = ones(size(new_branches, 1), 1);
                    for j=1:size(new_branches, 1)
                        time_enter = total_cost + new_branches(j, 7) - new_branches(j, 4) + abs(this_row - new_branches(j, 1)) + saver;
                        time_exit = total_cost + new_branches(j, 7) + saver;
                        if new_branches(j, 2) == 0
                            avoidance = [avoidance_map{new_branches(j, 1), 2:size_column-1}];
                        else
                            if new_branches(j, 3) == 1
                                avoidance = [avoidance_map{new_branches(j, 1), 2:new_branches(j, 2)}];
                            elseif new_branches(j, 3) == 2
                                avoidance = [avoidance_map{new_branches(j, 1), new_branches(j, 2):size_column-1}];
                            end
                        end
                        if ~isempty(fastintersect((time_enter:time_exit) + 1, avoidance + 1))
                            keeper(j) = 0;
                        end
                    end
                    new_branches = new_branches(logical(keeper), :);
                    %pick random branches
                    if randomized >= 1
                        if sum(new_branches(:, 9) > 0) > best_of
                            idxs = datasample([1:size(new_branches, 1)], best_of, 'Replace', false, 'Weights', new_branches(:, 9)./sum(new_branches(:, 9)));
                            new_branches = new_branches(idxs, :);
                        end
                    end
                    if size(new_branches, 1) > best_of
                        [vals, idxes] = maxk(new_branches(:, 9), best_of);
                        new_branches = new_branches(idxes, :);
                    end
                    %Save branches
                    if isempty(new_branches)
                        new_branches = branch(k, :);
                        new_branches(1, [4,5,6]) = 0;
                        new_branches(1, 10) = k;
                    end
                    horizon_branches = [horizon_branches; new_branches];
                    next_cumulative_reward_1{k} = temp_cumulative_reward_1;
                    next_cumulative_reward_2{k} = temp_cumulative_reward_2;
                end
                decision_tree{i} = horizon_branches;
            end
            %% summed heuristic
%             rewards = cell(horizon+deeper, 1);
%             costs = cell(horizon+deeper, 1);
%             leafs = decision_tree{end};
%             rewards{end} = leafs(:, 8);
%             costs{end} = leafs(:, 7);
%             for i=(horizon+deeper-1):-1:1
%                 prev_rewards = rewards{i+1};
%                 prev_costs = costs{i+1};
%                 temp_rewards = zeros(max(leafs(:, 10)), 1);
%                 temp_costs = temp_rewards;
%                 for j=1:size(leafs, 1)
%                     idx = leafs(j, 10);
%                     temp_rewards(idx) = temp_rewards(idx) + prev_rewards(j) + leafs(j, 8);
%                     temp_costs(idx) = temp_costs(idx) + prev_costs(j) + leafs(j, 7);
%                 end
%                 leafs = decision_tree{i};
%                 for j=1:size(leafs, 1)
%                     temp_rewards(j) = temp_rewards(j) + leafs(j, 8);
%                     temp_costs(j) = temp_costs(j) + leafs(j, 7);
%                 end
%                 rewards{i} = temp_rewards;
%                 costs{i} = temp_costs;
%             end
%             best_move = [];
%             prev_move = 1;
%             for i=1:horizon
%                 leafs = decision_tree{i};
%                 r = rewards{i};
%                 c = costs{i};
%                 number = 1:length(r);
%                 if ~isnan(sum(r)/sum(c)) && sum(r) > 0
%                     idxs = (leafs(:,10) == prev_move);
%                     leaves = leafs(idxs, :);
%                     r = r(idxs);
%                     c = c(idxs);
%                     number = number(idxs);
%                     [val, idx] = max(r./c, [], 'omitnan');
%                     best_move = [best_move; leaves(idx, :)];
%                     prev_move = number(idx);
%                 end
%             end
            %% find branch
            found = 0;
            best_move = [];
            for i=(horizon+deeper):-1:1
                leafs = decision_tree{i};
                if found
                    best_move = [leafs(best_move(1, 10), :); best_move];
                else
                    if ~isempty(leafs)
                        if any(leafs(:, 9) > 0)
%                             %random
%                             idx = randsample([1:size(leafs, 1)], 1, true, leafs(:, 9)./sum(leafs(:, 9)));
%                             best_heuristic = leafs(idx, 9);
                            %not random
                            [best_heuristic, idx] = max(leafs(:, 9));
                            if best_heuristic > 0
                                best_move = leafs(idx, :);
                                found = 1;
                            end
                        end
                    end
                end
            end
            %% make all moves
            if isempty(best_move)
                break;
            end
            best_move = best_move(1:min(num_moves, size(best_move, 1)), :);
            while ~isempty(best_move)
                if best_move(1, 4) == 0
                    break;
                end
                if best_move(1, 2) == 0
                    best_loop_heuristic = 1;
                    best_muted_heuristic = 0;
                    best_loop_index = best_move(1, 1);
                else
                    best_muted_heuristic = 1;
                    best_loop_heuristic = 0;
                    best_muted_row = best_move(1, 1);
                    best_muted_col = best_move(1, 2);
                end
                % add to tour
                if current_side == 1
                    if (best_loop_heuristic >= best_muted_heuristic)
                        if current_row <= best_loop_index
                            merp = 1;
                        else
                            merp = -1;
                        end
                        for i=current_row+merp:merp:best_loop_index
                                total_cost = total_cost + 1;
                                tour = [tour; (i-1)*size_column+1];
                                avoidance_map{i, 1} = [avoidance_map{i, 1}, total_cost+saver];
                                total_reward = total_reward + cumulative_reward_1(i, 1);
                                cumulative_reward_1(i, :) = max(0, cumulative_reward_1(i, :) - cumulative_reward_1(i, 1));
                                cumulative_reward_2(i, 1) = cumulative_reward_2(i, 2);
                        end
                        for i=2:size_column
                            total_cost = total_cost + 1;
                            tour = [tour; (best_loop_index-1)*size_column + i]; 
                            avoidance_map{best_loop_index, i} = [avoidance_map{best_loop_index, i}, total_cost+saver];
                        end
                        current_row = best_loop_index;
                        current_side = 2;
                        total_reward = total_reward + cumulative_reward_1(best_loop_index, size_column);
                        cumulative_reward_1(best_loop_index, :) = 0;
                        cumulative_reward_2(best_loop_index, :) = 0;
                    else
                        if current_row <= best_muted_row
                            merp = 1;
                        else
                            merp = -1;
                        end
                        for i=current_row+merp:merp:best_muted_row
                                total_cost = total_cost + 1;
                                tour = [tour; (i-1)*size_column+1];
                                avoidance_map{i, 1} = [avoidance_map{i, 1}, total_cost+saver];
                                total_reward = total_reward + cumulative_reward_1(i, 1);
                                cumulative_reward_1(i, :) = max(0, cumulative_reward_1(i, :) - cumulative_reward_1(i, 1));
                                cumulative_reward_2(i, 1) = cumulative_reward_2(i, 2);
                        end
                        for i=2:best_muted_col
                            total_cost = total_cost + 1;
                            tour = [tour; (best_muted_row-1)*size_column + i];
                            avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, total_cost+saver];
                        end
                        for i=best_muted_col-1:-1:1
                            total_cost = total_cost + 1;
                            tour = [tour; (best_muted_row-1)*size_column + i];
                            avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, total_cost+saver];
                        end
                        accumulated_reward = cumulative_reward_1(best_muted_row, best_muted_col);
                        total_reward = total_reward + accumulated_reward;
                        current_row = best_muted_row;
                        cumulative_reward_1(best_muted_row, :) = max(0, cumulative_reward_1(best_muted_row, :) - accumulated_reward);
                        cumulative_reward_2(best_muted_row, 1:best_muted_col) = cumulative_reward_2(best_muted_row, min(best_muted_col + 1, size_column));
                        if best_muted_col == size_column
                            temp_cumulative_reward_2(best_muted_row, :) = 0;
                        end
                    end
                elseif current_side == 2
                    if (best_loop_heuristic >= best_muted_heuristic)
                        if current_row <= best_loop_index
                            merp = 1;
                        else
                            merp = -1;
                        end
                        for i=current_row+merp:merp:best_loop_index
                                total_cost = total_cost + 1;
                                tour = [tour; (i)*size_column];
                                avoidance_map{i, size_column} = [avoidance_map{i, size_column}, total_cost+saver];
                                total_reward = total_reward + cumulative_reward_2(i, size_column);
                                cumulative_reward_2(i, :) = max(0, cumulative_reward_2(i, :) - cumulative_reward_2(i, size_column));
                                cumulative_reward_1(i, size_column) = cumulative_reward_1(i, size_column-1);
                        end
                        for i=size_column-1:-1:1
                            total_cost = total_cost + 1;
                            tour = [tour; (best_loop_index-1)*size_column + i]; 
                            avoidance_map{best_loop_index, i} = [avoidance_map{best_loop_index, i}, total_cost+saver];
                        end
                        current_row = best_loop_index;
                        current_side = 1;
                        total_reward = total_reward + cumulative_reward_2(best_loop_index, 1);
                        cumulative_reward_2(best_loop_index, :) = 0;
                        cumulative_reward_1(best_loop_index, :) = 0;
                    else
                        if current_row <= best_muted_row
                            merp = 1;
                        else
                            merp = -1;
                        end
                        for i=current_row+merp:merp:best_muted_row
                                total_cost = total_cost + 1;
                                tour = [tour; (i)*size_column];
                                avoidance_map{i, size_column} = [avoidance_map{i, size_column}, total_cost+saver];
                                total_reward = total_reward + cumulative_reward_2(i, size_column);
                                cumulative_reward_2(i, :) = max(0, cumulative_reward_2(i, :) - cumulative_reward_2(i, size_column));
                                cumulative_reward_1(i, size_column) = cumulative_reward_1(i, size_column-1);
                        end
                        for i=size_column-1:-1:best_muted_col
                            total_cost = total_cost + 1;
                            tour = [tour; (best_muted_row-1)*size_column + i];
                            avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, total_cost+saver];
                        end
                        for i=best_muted_col+1:1:size_column
                            total_cost = total_cost + 1;
                            tour = [tour; (best_muted_row-1)*size_column + i];
                            avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, total_cost+saver];
                        end
                        accumulated_reward = cumulative_reward_2(best_muted_row, best_muted_col);
                        total_reward = total_reward + accumulated_reward;
                        current_row = best_muted_row;
                        cumulative_reward_2(best_muted_row, :) = max(0, cumulative_reward_2(best_muted_row, :) - accumulated_reward);
                        cumulative_reward_1(best_muted_row, size_column:-1:best_muted_col) = cumulative_reward_1(best_muted_row, max(best_muted_col - 1, 1));
                        if best_muted_col == size_column
                            temp_cumulative_reward_1(best_muted_row, :) = 0;
                        end
                    end
                end
                best_move(1, :) = [];
            end
        end
        %% Finish tour
        if current_side == 1
            if current_row > ending_row
                derp = -1;
            else
                derp = 1;
            end
            for i=current_row+derp:derp:ending_row
                total_cost = total_cost + 1;
                tour = [tour; (i-1)*size_column+1];
                avoidance_map{i, 1} = [avoidance_map{i, 1}, total_cost+saver];
                total_reward = total_reward + cumulative_reward_1(i, 1);
                cumulative_reward_1(i, :) = max(0, cumulative_reward_1(i, :) - cumulative_reward_1(i, 1));
                cumulative_reward_2(i, 1) = cumulative_reward_2(i, 2);
            end
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
                home_stretch = row_distance*abs((1:size_row)' - current_row) + cumulative_cost_2(:,1) + row_distance*abs((1:size_row)' - ending_row);
                for l=1:size_row
                    [branch, k] = min(home_stretch);
                    row_avoidance = [avoidance_map{k, 2:size_column-1}];
                    to_loop_cost = row_distance*abs(current_row - k);
                    through_loop_cost = cumulative_cost_2(k, 1);
                    to_end_cost = row_distance*abs(k - ending_row);
                    time_enter = total_cost + to_loop_cost + saver;
                    time_exit = total_cost + to_loop_cost + through_loop_cost + saver;
                    if isempty(fastintersect(time_enter:time_exit, row_avoidance))
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
                        avoidance_map{i, size_column} = [avoidance_map{i, size_column}, total_cost+saver];
                        total_reward = total_reward + cumulative_reward_2(i, size_column);
                        cumulative_reward_2(i, :) = max(0, cumulative_reward_2(i, :) - cumulative_reward_2(i, size_column));
                        cumulative_reward_1(i, size_column) = cumulative_reward_1(i, size_column-1);
                    end
                    for i=size_column-1:-1:1
                        total_cost = total_cost + 1;
                        tour = [tour; (k-1)*size_column + i];
                        avoidance_map{k, i} = [avoidance_map{k, i}, total_cost+saver];
                    end
                    total_reward = total_reward + cumulative_reward_2(k, 1);
                    if current_row > ending_row
                        derp = -1;
                    else
                        derp = 1;
                    end
                    for i=k+derp:derp:ending_row
                        total_cost = total_cost + 1;
                        tour = [tour; (i-1)*size_column+1];
                        avoidance_map{i, 1} = [avoidance_map{i, 1}, total_cost+saver];
                        total_reward = total_reward + cumulative_reward_1(i, 1);
                        cumulative_reward_1(i, :) = max(0, cumulative_reward_1(i, :) - cumulative_reward_1(i, 1));
                        cumulative_reward_2(i, 1) = cumulative_reward_2(i, 2);
                    end
                    current_row = k;
                else
                    total_cost = total_cost + 1;
                    waiting = waiting - 1;
                    tour = [tour; tour(end)];
                end
            end
        end
        if total_reward > sum(sum(old_reward_map))
            %pause
        end
        %% Check if satisfied
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
end


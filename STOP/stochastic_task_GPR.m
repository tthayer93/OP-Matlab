function [ energy_cost, resource_cost, collected_reward, wasted_resources, tour, priority_map, avoidance_map ] = stochastic_task_GPR( row_distance, vine_distance, priority_map, service_cost_map, beginning, ending, energy_budget, resource_budget, avoidance_map, saver )
%STOCHASTIC_TASK_GPR Adaption of GPR for simulating vertex service cost
%
%	Version: 2.0
%	Date: 15/12/20
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function is an adaption of the Greedy Partial Row (GPR) heuristic to simulate visiting vertices with unknown service costs, according to the problem designated in "Task Planning on Stochastic Aisle Graphs for Precision Agriculture" by Kan, Thayer, Carpin, and Karydis.
%	Inputs: 
%		row_distance: Constant distance between each vine row
%		vine_distance: Constant distance between each vine in a row
%		priority_map: A matrix containing the priority multipliers for every vertex (analogous to the reward_map in the original GPR heuristic)
%		service_cost_map: A matrix containing the actual service costs for every vertex, which are unknown to the robot and therefore the route planner
%		beginning: The vertex at which the route should begin (should be on left outside column)
%		ending: The vertex at which the route should end (should be on left outside column)
%		energy_budget: The budget defining the total allowable movement cost for the robot (analogous to the budget in the original GPR heuristic)
%		resource_budget: The budget defining the total allowable vertex service cost for the robot
%		avoidance_map: For avoiding other robots in the team scenario
%		saver: Initial value for saver, use if extending a previously built tour
%	Outputs:
%		energy_cost: The total movement cost for the computed route
%		resource_cost: The total vertex service cost for the computed route
%		collected_reward: The total vertex service reward for each vertex in the computed route, where the reward for each vertex is its priority multiplier times its service cost. This vector should be the same length as the tour with a 1 to 1 correspondance
%		wasted_resources: The number of resources that were wasted when trying to service a vertex unsuccessfully
%		tour: The computed route
%		priority_map: An updated priority map after servicing vertices in the computed route
%		avoidance_map: For avoiding other robots in the team scenario
	
	row_distance = 1;
	vine_distance = 1;
	if(nargin < 10)
		saver = 0;
	end
    if(nargin < 9)
        avoidance_map = cell(size(priority_map));
    end
    initial_reward_map = priority_map;
    avoidance_map_initial = avoidance_map;
    unsatisfied = 1;
    while unsatisfied
        %% Setup initial values
		waiting = saver;
		collected_reward = [];
		wasted_resources = [];
		priority_map = initial_reward_map;
        avoidance_map = avoidance_map_initial;
        energy_cost = 0;
		resource_cost = 0;
		beginning_row = ceil(beginning/size(priority_map,2));
        ending_row = ceil(ending/size(priority_map,2));
		current_row = beginning_row;
        current_side = 1;
        tour = [beginning*ones(waiting+1, 1)];
        collected_reward = [priority_map(current_row, 1) * service_cost_map(current_row, 1); zeros(waiting, 1)];
		wasted_resources = [zeros(waiting+1, 1)];
		resource_cost = resource_cost + double(logical(priority_map(current_row, 1))) * service_cost_map(current_row, 1);
        priority_map(current_row, 1) = 0;
        cumulative_reward_1 = zeros(size(priority_map));
        cumulative_reward_2 = zeros(size(priority_map));
        cumulative_cost_1 = zeros(size(priority_map));
        cumulative_cost_2 = zeros(size(priority_map));
        size_row = size(priority_map, 1); 
        size_column = size(priority_map, 2);
        for i=1:size_row
            for j=1:size_column
                cumulative_reward_1(i, j) = sum(priority_map(i, 1:j));
                cumulative_reward_2(i, size_column - j + 1) = sum(priority_map(i, size_column:-1:size_column-j+1));
                cumulative_cost_1(i, j) = (j - 1) * vine_distance;
                cumulative_cost_2(i, size_column - j + 1) = (j - 1) * vine_distance;
            end
        end

        %% Find path
        while (sum(sum(cumulative_reward_1)) > 0 || sum(sum(cumulative_reward_2)) > 0) && (resource_budget - resource_cost > 0)
            if current_side == 1
				%% If on left side of vineyard
                distance_side = row_distance * abs(current_row - (1:size_row)');
                loop_heuristic = cumulative_reward_1(:, size_column)./(cumulative_cost_1(:, size_column) + distance_side);
                muted_loop_heuristic = cumulative_reward_1./(2 * cumulative_cost_1 + distance_side);
				%% Find best row
                unsatisfied1 = 0;
                satisfied = 0;
                idk = 0;
                while ~satisfied
                    [best_loop_heuristic, best_loop_index] = max(loop_heuristic);
                    row_avoidance = [avoidance_map{best_loop_index, 2:size_column-1}];
                    to_loop_cost = row_distance*abs(current_row - best_loop_index);
                    through_loop_cost = cumulative_cost_1(best_loop_index, size_column);
                    to_end_cost = row_distance*abs(best_loop_index - ending_row);
                    time_enter = energy_cost + to_loop_cost + waiting;
                    time_exit = time_enter + through_loop_cost + waiting;
                    if any(ismember(time_enter:time_exit, row_avoidance))
                        if (loop_heuristic(best_loop_index) > 0) && (to_loop_cost + 2*through_loop_cost + to_end_cost + 1 <= energy_budget - energy_cost)
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
                %% Find best partial row
				unsatisfied2 = 0;
                satisfied = 0;
                while ~satisfied
                    best_muted_heuristic = max(muted_loop_heuristic(:));
                    [best_muted_row, best_muted_col] = find(muted_loop_heuristic == best_muted_heuristic, 1);
                    muted_avoidance = [avoidance_map{best_muted_row, 2:best_muted_col}];
                    to_muted_cost = row_distance*abs(current_row - best_muted_row);
                    through_muted_cost = 2*cumulative_cost_1(best_muted_row, best_muted_col);
                    to_end_cost = row_distance*abs(best_muted_row - ending_row);
                    time_enter = energy_cost + to_muted_cost + waiting;
                    time_exit = time_enter + through_muted_cost + waiting;
                    if any(ismember(time_enter:time_exit, muted_avoidance))
                        if (muted_loop_heuristic(best_muted_row, best_muted_col) > 0) && (to_muted_cost + through_muted_cost + to_end_cost + 1 <= energy_budget - energy_cost)
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
                %% If no available movements then wait
				if unsatisfied1 && unsatisfied2
                    if idk > 0
                        %energy_cost = energy_cost + 1;
						waiting = waiting + 1;
                        tour = [tour; tour(end)];
						collected_reward = [collected_reward; 0];
						wasted_resources = [wasted_resources; 0];
						continue;
                    else
                        break;
                    end
				end
				%% idk
                %[best_loop_heuristic, best_loop_index] = max(loop_heuristic);
                %best_muted_heuristic = max(muted_loop_heuristic(:));
                %[best_muted_row, best_muted_col] = find(muted_loop_heuristic == best_muted_heuristic, 1);
                %loop_cost = row_distance*abs(current_row - best_loop_index) + 2*cumulative_cost_1(best_loop_index, size_column) + row_distance*abs(best_loop_index - ending_row);
                %muted_loop_cost = row_distance*abs(current_row - best_muted_row) + 2*cumulative_cost_1(best_muted_row, best_muted_col) + row_distance*abs(best_muted_row - ending_row);
				if (best_loop_heuristic >= best_muted_heuristic) && (loop_cost <= (energy_budget - energy_cost))
					%% Use best row
                    %total_cost = total_cost + row_distance*abs(current_row - best_loop_index) + cumulative_cost_1(best_loop_index, size_column);
                    %total_reward = total_reward + cumulative_reward_1(best_loop_index, size_column);
                    %tour = [tour; best_loop_index, size_column];
                    if current_row <= best_loop_index
                        merp = 1;
                    else
                        merp = -1;
					end
                    for i=current_row+merp:merp:best_loop_index
                            energy_cost = energy_cost + 1;
                            tour = [tour; (i-1)*size_column+1];
                            avoidance_map{i, 1} = [avoidance_map{i, 1}, energy_cost + waiting];
							resource_cost = resource_cost + double(logical(priority_map(i, 1))) * service_cost_map(i, 1);
							if resource_cost > resource_budget
								collected_reward = [collected_reward; 0];
								wasted_resources = [wasted_resources; resource_budget - (resource_cost - double(logical(priority_map(i, 1))) * service_cost_map(i, 1))];
								break;
							end
							collected_reward = [collected_reward; priority_map(i, 1) * service_cost_map(i, 1)];
							wasted_resources = [wasted_resources; 0];
							priority_map(i, 1) = 0;
                            %total_reward = total_reward + cumulative_reward_1(i, 1);
							cumulative_reward_1(i, :) = max(0, cumulative_reward_1(i, :) - cumulative_reward_1(i, 1));
                            cumulative_reward_2(i, 1) = cumulative_reward_2(i, 2);
					end
					current_row = best_loop_index;
					if resource_cost >= resource_budget
						break;
					end
                    for i=2:size_column
                        energy_cost = energy_cost + 1;
                        tour = [tour; (best_loop_index-1)*size_column + i]; 
                        avoidance_map{best_loop_index, i} = [avoidance_map{best_loop_index, i}, energy_cost + waiting];
						resource_cost = resource_cost + double(logical(priority_map(current_row, i))) * service_cost_map(current_row, i);
						if resource_cost <= resource_budget
							collected_reward = [collected_reward; priority_map(current_row, i) * service_cost_map(current_row, i)];
							wasted_resources = [wasted_resources; 0];
							priority_map(current_row, i) = 0;
						else
							collected_reward = [collected_reward; 0];
							if sum(wasted_resources) == 0
								wasted_resources = [wasted_resources; resource_budget - (resource_cost - double(logical(priority_map(current_row, i))) * service_cost_map(current_row, i))];
							else
								wasted_resources = [wasted_resources; 0];
							end
							resource_cost = resource_budget;
						end
					end
                    current_side = 2;
                    %total_reward = total_reward + cumulative_reward_1(best_loop_index, size_column);
                    cumulative_reward_1(best_loop_index, :) = 0;
                    cumulative_reward_2(best_loop_index, :) = 0;
				elseif (muted_loop_cost <= (energy_budget - energy_cost))
					%% Use best partial row
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
                            energy_cost = energy_cost + 1;
                            tour = [tour; (i-1)*size_column+1];
                            avoidance_map{i, 1} = [avoidance_map{i, 1}, energy_cost + waiting];
							resource_cost = resource_cost + double(logical(priority_map(i, 1))) * service_cost_map(i, 1);
							if resource_cost > resource_budget
								collected_reward = [collected_reward; 0];
								wasted_resources = [wasted_resources; resource_budget - (resource_cost - double(logical(priority_map(i, 1))) * service_cost_map(i, 1))];
								break;
							end
							collected_reward = [collected_reward; priority_map(i, 1) * service_cost_map(i, 1)];
							wasted_resources = [wasted_resources; 0];
							priority_map(i, 1) = 0;
							%total_reward = total_reward + cumulative_reward_1(i, 1);
                            cumulative_reward_1(i, :) = max(0, cumulative_reward_1(i, :) - cumulative_reward_1(i, 1));
                            cumulative_reward_2(i, 1) = cumulative_reward_2(i, 2);
					end
					current_row = best_muted_row;
					if resource_cost >= resource_budget
						break;
					end
                    for i=2:best_muted_col
                        energy_cost = energy_cost + 1;
                        tour = [tour; (best_muted_row-1)*size_column + i];
                        avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, energy_cost + waiting];
						resource_cost = resource_cost + double(logical(priority_map(current_row, i))) * service_cost_map(current_row, i);
						if resource_cost <= resource_budget
							collected_reward = [collected_reward; priority_map(current_row, i) * service_cost_map(current_row, i)];
							wasted_resources = [wasted_resources; 0];
							priority_map(current_row, i) = 0;
						else
							collected_reward = [collected_reward; 0];
							if sum(wasted_resources) == 0
								wasted_resources = [wasted_resources; resource_budget - (resource_cost - double(logical(priority_map(current_row, i))) * service_cost_map(current_row, i))];
							else
								wasted_resources = [wasted_resources; 0];
							end
							resource_cost = resource_budget;
						end
                    end
                    for i=best_muted_col-1:-1:1
                        energy_cost = energy_cost + 1;
                        tour = [tour; (best_muted_row-1)*size_column + i];
						collected_reward = [collected_reward; 0];
						wasted_resources = [wasted_resources; 0];
                        avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, energy_cost + waiting];
                    end
                    accumulated_reward = cumulative_reward_1(best_muted_row, best_muted_col);
                    %total_reward = total_reward + accumulated_reward;
                    cumulative_reward_1(best_muted_row, :) = max(0, cumulative_reward_1(best_muted_row, :) - accumulated_reward);
                    %cumulative_reward_2(best_muted_row, 1:best_muted_col) = cumulative_reward_2(best_muted_row, best_muted_col + 1);
                    cumulative_reward_2(best_muted_row, 1:best_muted_col) = cumulative_reward_2(best_muted_row, min(best_muted_col + 1, size_column));
				else
					%% Best row and partial row are infeasible
                    temp_cost = row_distance*abs(current_row - (1:size_row))' + 2*cumulative_cost_1 + row_distance*abs((1:size_row) - ending_row)';
                    reachable = (temp_cost <= (energy_budget - energy_cost));
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
				%% If on right side of vineyard
                distance_side = row_distance * abs(current_row - (1:size_row)');
                loop_heuristic = cumulative_reward_2(:, 1)./(cumulative_cost_2(:, 1) + distance_side);
                muted_loop_heuristic = cumulative_reward_2./(2 * cumulative_cost_2 + distance_side);
				%% Find best row
                unsatisfied1 = 0;
                satisfied = 0;
                idk = 0;
                while ~satisfied
                    [best_loop_heuristic, best_loop_index] = max(loop_heuristic);
                    row_avoidance = [avoidance_map{best_loop_index, 2:size_column-1}];
                    to_loop_cost = row_distance*abs(current_row - best_loop_index);
                    through_loop_cost = cumulative_cost_2(best_loop_index, 1);
                    to_end_cost = row_distance*abs(best_loop_index - ending_row);
                    time_enter = energy_cost + to_loop_cost + waiting;
                    time_exit = time_enter + through_loop_cost + waiting;
                    if any(ismember(time_enter:time_exit, row_avoidance))
                        if (loop_heuristic(best_loop_index) > 0) && (to_loop_cost + through_loop_cost + to_end_cost + 1 <= energy_budget - energy_cost)
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
                %% Find best partial row
				unsatisfied2 = 0;
                satisfied = 0;
                while ~satisfied
                    best_muted_heuristic = max(muted_loop_heuristic(:));
                    [best_muted_row, best_muted_col] = find(muted_loop_heuristic == best_muted_heuristic, 1);
                    muted_avoidance = [avoidance_map{best_muted_row, best_muted_col:size_column-1}];
                    to_muted_cost = row_distance*abs(current_row - best_muted_row);
                    through_muted_cost = 2*cumulative_cost_2(best_muted_row, best_muted_col);
                    to_end_cost = row_distance*abs(best_muted_row - ending_row) + cumulative_cost_2(ending_row, 1);
                    time_enter = energy_cost + to_muted_cost + waiting;
                    time_exit = time_enter + through_muted_cost + waiting;
                    if any(ismember(time_enter:time_exit, muted_avoidance))
                        if (muted_loop_heuristic(best_muted_row, best_muted_col) > 0) && (to_muted_cost + through_muted_cost + to_end_cost + 1 <= energy_budget - energy_cost)
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
                %% If no available movements then wait
				if unsatisfied1 && unsatisfied2
                    if idk > 0
                        %energy_cost = energy_cost + 1;
						waiting = waiting + 1;
                        tour = [tour; tour(end)];
						collected_reward = [collected_reward; 0];
						wasted_resources = [wasted_resources; 0];
						continue;
                    else
                        break;
                    end
				end
				%% idk
                %[best_loop_heuristic, best_loop_index] = max(loop_heuristic);
                %best_muted_heuristic = max(muted_loop_heuristic(:));
                %[best_muted_row, best_muted_col] = find(muted_loop_heuristic == best_muted_heuristic, 1);
                %loop_cost = row_distance*abs(current_row - best_loop_index) + cumulative_cost_2(best_loop_index, 1) + row_distance*abs(best_loop_index - ending_row);
                %muted_loop_cost = row_distance*abs(current_row - best_muted_row) + 2*cumulative_cost_2(best_muted_row, best_muted_col) + row_distance*abs(best_muted_row - ending_row) + cumulative_cost_2(ending_row, 1);
                if (best_loop_heuristic >= best_muted_heuristic) && (loop_cost <= (energy_budget - energy_cost))
					%% Use best row
                    %total_cost = total_cost + row_distance*abs(current_row - best_loop_index) + cumulative_cost_2(best_loop_index, 1);
                    %total_reward = total_reward + cumulative_reward_2(best_loop_index, 1);
                    %tour = [tour; best_loop_index, 1];
                    if current_row <= best_loop_index
                        merp = 1;
                    else
                        merp = -1;
                    end
                    for i=current_row+merp:merp:best_loop_index
                            energy_cost = energy_cost + 1;
                            tour = [tour; (i)*size_column];
                            avoidance_map{i, size_column} = [avoidance_map{i, size_column}, energy_cost + waiting];
                            resource_cost = resource_cost + double(logical(priority_map(i, size_column))) * service_cost_map(i, size_column);
							if resource_cost > resource_budget
								collected_reward = [collected_reward; 0];
								wasted_resources = [wasted_resouces; resource_budget - (resource_cost - double(logical(priority_map(i, size_column))) * service_cost_map(i, size_column))];
								break;
							end
							collected_reward = [collected_reward; priority_map(i, size_column) * service_cost_map(i, size_column)];
							wasted_resources = [wasted_resources; 0];
							priority_map(i, size_column) = 0;
							%total_reward = total_reward + cumulative_reward_2(i, size_column);
                            cumulative_reward_2(i, :) = max(0, cumulative_reward_2(i, :) - cumulative_reward_2(i, size_column));
                            cumulative_reward_1(i, size_column) = cumulative_reward_1(i, size_column-1);
					end
					current_row = best_loop_index;
					if resource_cost >= resource_budget
						break;
					end
                    for i=size_column-1:-1:1
                        energy_cost = energy_cost + 1;
                        tour = [tour; (best_loop_index-1)*size_column + i]; 
                        avoidance_map{best_loop_index, i} = [avoidance_map{best_loop_index, i}, energy_cost + waiting];
						resource_cost = resource_cost + double(logical(priority_map(current_row, i))) * service_cost_map(current_row, i);
						if resource_cost <= resource_budget
							collected_reward = [collected_reward; priority_map(current_row, i) * service_cost_map(current_row, i)];
							wasted_resources = [wasted_resources; 0];
							priority_map(current_row, i) = 0;
						else
							collected_reward = [collected_reward; 0];
							if sum(wasted_resources) == 0
								wasted_resources = [wasted_resources; resource_budget - (resource_cost - double(logical(priority_map(current_row, i))) * service_cost_map(current_row, i))];
							else
								wasted_resources = [wasted_resources; 0];
							end
							resource_cost = resource_budget;
						end
					end
                    current_side = 1;
                    %total_reward = total_reward + cumulative_reward_2(best_loop_index, 1);
                    cumulative_reward_2(best_loop_index, :) = 0;
                    cumulative_reward_1(best_loop_index, :) = 0;
                elseif (muted_loop_cost <= (energy_budget - energy_cost))
					%% Use best partial row
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
                            energy_cost = energy_cost + 1;
                            tour = [tour; (i)*size_column];
                            avoidance_map{i, size_column} = [avoidance_map{i, size_column}, energy_cost + waiting];
                            resource_cost = resource_cost + double(logical(priority_map(i, size_column))) * service_cost_map(i, size_column);
							if resource_cost > resource_budget
								collected_reward = [collected_reward; 0];
								wasted_resources = [wasted_resources; resource_budget - (resource_cost - double(logical(priority_map(i, size_column))) * service_cost_map(i, size_column))];
								break;
							end
							collected_reward = [collected_reward; priority_map(i, size_column) * service_cost_map(i, size_column)];
							wasted_resources = [wasted_resources; 0];
							priority_map(i, size_column) = 0;
							%total_reward = total_reward + cumulative_reward_2(i, size_column);
                            cumulative_reward_2(i, :) = max(0, cumulative_reward_2(i, :) - cumulative_reward_2(i, size_column));
                            cumulative_reward_1(i, size_column) = cumulative_reward_1(i, size_column-1);
					end
					current_row = best_muted_row;
					if resource_cost >= resource_budget
						break;
					end
                    for i=size_column-1:-1:best_muted_col
                        energy_cost = energy_cost + 1;
                        tour = [tour; (best_muted_row-1)*size_column + i];
                        avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, energy_cost + waiting];
						resource_cost = resource_cost + double(logical(priority_map(current_row, i))) * service_cost_map(current_row, i);
						if resource_cost <= resource_budget
							collected_reward = [collected_reward; priority_map(current_row, i) * service_cost_map(current_row, i)];
							wasted_resources = [wasted_resources; 0];
							priority_map(current_row, i) = 0;
						else
							collected_reward = [collected_reward; 0];
							if sum(wasted_resources) == 0
								wasted_resources = [wasted_resources; resource_budget - (resource_cost - double(logical(priority_map(current_row, i))) * service_cost_map(current_row, i))];
							else
								wasted_resources = [wasted_resources; 0];
							end
							resource_cost = resource_budget;
						end
                    end
                    for i=best_muted_col+1:1:size_column
                        energy_cost = energy_cost + 1;
                        tour = [tour; (best_muted_row-1)*size_column + i];
						collected_reward = [collected_reward; 0];
						wasted_resources = [wasted_resources; 0];
                        avoidance_map{best_muted_row, i} = [avoidance_map{best_muted_row, i}, energy_cost + waiting];
                    end
                    accumulated_reward = cumulative_reward_2(best_muted_row, best_muted_col);
                    %total_reward = total_reward + accumulated_reward;
                    cumulative_reward_2(best_muted_row, :) = max(0, cumulative_reward_2(best_muted_row, :) - accumulated_reward);
                    %cumulative_reward_1(best_muted_row, size_column:-1:best_muted_col) = cumulative_reward_1(best_muted_row, best_muted_col - 1);
                    cumulative_reward_1(best_muted_row, size_column:-1:best_muted_col) = cumulative_reward_1(best_muted_row, max(best_muted_col - 1, 1));
				else
					%% Best row and partial row are infeasible
                    temp_cost = row_distance*abs(current_row - (1:size_row))' + cumulative_cost_2(:, 1) + row_distance*abs((1:size_row) - ending_row)';
                    reachable = (temp_cost <= (energy_budget - energy_cost));
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
		%% Go to end vertex
        if current_side == 1
            %total_cost = total_cost + row_distance*abs(current_row - ending_row);
            %tour = [tour; ending_row, 1];
            if current_row > ending_row
                derp = -1;
            else
                derp = 1;
            end
            for i=current_row+derp:derp:ending_row
                energy_cost = energy_cost + 1;
                tour = [tour; (i-1)*size_column+1];
                avoidance_map{i, 1} = [avoidance_map{i, 1}, energy_cost + waiting];
                %total_reward = total_reward + cumulative_reward_1(i,1);
				resource_cost = resource_cost + double(logical(priority_map(i, 1))) * service_cost_map(i, 1);
				if resource_cost > resource_budget
					collected_reward = [collected_reward; 0];
					if sum(wasted_resources) == 0
						wasted_resources = [wasted_resources; resource_budget - (resource_cost - double(logical(priority_map(i, 1))) * service_cost_map(i, 1))];
					else
						wasted_resources = [wasted_resources; 0];
					end
					resource_cost = resource_budget;
				else
					collected_reward = [collected_reward; priority_map(i, 1) * service_cost_map(i, 1)];
					wasted_resources = [wasted_resources; 0];
					priority_map(i, 1) = 0;
				end
            end
            %tour = [tour; ending];
            current_row = ending_row;
        end
        unsatisfied = 0;
        if current_side == 2 % Will only be at side 2 if the row needed to reach the ending has already been traversed
            %total_cost = total_cost + row_distance*abs(current_row - ending_row) + cumulative_cost_2(current_row, 1);
            %tour = [tour; ending_row, size_column; ending_row, 1];

            while tour(end) ~= ending
                if energy_budget - energy_cost < cumulative_cost_2(1,1)
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
                    time_enter = energy_cost + to_loop_cost + waiting;
                    time_exit = time_enter + through_loop_cost + waiting;
                    if ~any(ismember(time_enter:time_exit, row_avoidance))
                        satisfied = 1;
                        break;
                    else
                        home_stretch(k) = inf;
                    end
                end
                loop_cost = to_loop_cost + through_loop_cost + to_end_cost;
                if loop_cost <= (energy_budget - energy_cost) && satisfied
                    if current_row > k
                        derp = -1;
                    else
                        derp = 1;
                    end
                    for i=current_row+derp:derp:k
                        energy_cost = energy_cost + 1;
                        tour = [tour; (i)*size_column];
                        avoidance_map{i, size_column} = [avoidance_map{i, size_column}, energy_cost + waiting];
                        %total_reward = total_reward + cumulative_reward_1(i,1);
						resource_cost = resource_cost + double(logical(priority_map(i, size_column))) * service_cost_map(i, size_column);
						if resource_cost > resource_budget
							collected_reward = [collected_reward; 0];
							if sum(wasted_resources) == 0
								wasted_resources = [wasted_resources; resource_budget - (resource_cost - double(logical(priority_map(i, size_column))) * service_cost_map(i, size_column))];
							else
								wasted_resources = [wasted_resources; 0];
							end
							resource_cost = resource_budget;
						else
							collected_reward = [collected_reward; priority_map(i, size_column) * service_cost_map(i, size_column)];
							wasted_resources = [wasted_resources; 0];
							priority_map(i, size_column) = 0;
						end
                    end
                    for i=size_column-1:-1:1
                        energy_cost = energy_cost + 1;
                        tour = [tour; (k-1)*size_column + i];
                        avoidance_map{k, i} = [avoidance_map{k, i}, energy_cost + waiting];
						resource_cost = resource_cost + double(logical(priority_map(k, i))) * service_cost_map(k, i);
						if resource_cost > resource_budget
							collected_reward = [collected_reward; 0];
							if sum(wasted_resources) == 0
								wasted_resources = [wasted_resources; resource_budget - (resource_cost - double(logical(priority_map(k, i))) * service_cost_map(k, i))];
							else
								wasted_resources = [wasted_resources; 0];
							end
							resource_cost = resource_budget;
						else
							collected_reward = [collected_reward; priority_map(k, i) * service_cost_map(k, i)];
							wasted_resources = [wasted_resources; 0];
							priority_map(k, i) = 0;
						end
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
                        energy_cost = energy_cost + 1;
                        tour = [tour; (i-1)*size_column+1];
                        avoidance_map{i, 1} = [avoidance_map{i, 1}, energy_cost + waiting];
                        %total_reward = total_reward + cumulative_reward_1(i,1);
						resource_cost = resource_cost + double(logical(priority_map(i, 1))) * service_cost_map(i, 1);
						if resource_cost > resource_budget
							collected_reward = [collected_reward; 0];
							if sum(wasted_resources) == 0
								wasted_resources = [wasted_resouces; resource_budget - (resource_cost - double(logical(priority_map(i, 1))) * service_cost_map(i, 1))];
							else
								wasted_resources = [wasted_resources; 0];
							end
							resource_cost = resource_budget;
						else
							collected_reward = [collected_reward; priority_map(i, 1) * service_cost_map(i, 1)];
							wasted_resources = [wasted_resources; 0];
							priority_map(i, 1) = 0;
						end
                    end
                    current_row = k;
                    %tour = [tour; ending_row*size_column; ending];
                else
                    energy_cost = energy_cost + 1;
                    tour = [tour; tour(end)];
					collected_reward = [collected_reward; 0];
					wasted_resources = [wasted_resources; 0];
                end
            end
        end
        if ~unsatisfied
            break;
        end
        saver = saver + 1;
	end
    
end
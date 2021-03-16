function [ total_energy_costs, total_resource_costs, total_collected_rewards, total_wasted_resources, total_tours ] = stochastic_task_team_GPR_series( row_distance, vine_distance, priority_map, service_cost_map, beginning, ending, energy_budget, resource_budget, num_robots )
%STOCHASTIC_TASK_TEAM_GPR_SERIES Adaption of series GPR for simulating vertex service cost with a team of robots
%
%	Version: 1.0
%	Date: 17/12/20
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function is an adaption of the series Greedy Partial Row (sGPR) heuristic to simulate a team of robots visiting vertices with unknown service costs, according to the problem designated in "Task Planning on Stochastic Aisle Graphs for Precision Agriculture" by Kan, Thayer, Carpin, and Karydis. Routing continues until all vertices have been serviced, preferentially utilizing the robot with the shortest overall tour when extending a tour.
%	Inputs:
%		row_distance: Constant distance between each vine row
%		vine_distance: Constant distance between each vine in a row
%		priority_map: A matrix containing the priority multipliers for every vertex (analogous to the reward_map in the original GPR heuristic)
%		service_cost_map: A matrix containing the actual service costs for every vertex, which are unknown to the robot and therefore the route planner
%		beginning: The vertex at which the route should begin (should be on left outside column)
%		ending: The vertex at which the route should end (should be on left outside column)
%		energy_budget: The budget defining the total allowable movement cost for each robot (analogous to the budget in the original GPR heuristic)
%		resource_budget: The budget defining the total allowable vertex service cost for each robot
%		avoidance_map: A cell matrix the same size as priority_map, where each cell contains a list of the times during which a vertex is occupied
%		num_robots: The number of robots working together
%	Outputs:
%		total_energy_costs: A vector containing the total movement cost for each computed route
%		total_resource_costs: A vector containing the total vertex service cost for each computed route
%		total_collected_rewards: A cell vector where each cell contains the total vertex service reward for every vertex in the computed route of a single robot, where the reward for each vertex is its priority multiplier times its service cost. Each vector should be the same length as the tour with a 1 to 1 correspondance
%		total_wasted_resources: A cell vector where each cell contains the number of resources that were wasted when trying to service a vertex unsuccessfully for the route of a single robot
%		total_tours: A cell vector where each cell contains the computed route for a robot

    total_energy_costs = zeros(1, num_robots);
	total_resource_costs = zeros(1, num_robots);
    total_collected_rewards = cell(1, num_robots);
	total_wasted_resources = cell(1, num_robots);
    total_tours = cell(1, num_robots);
    avoidance_map = cell(size(priority_map));
	while sum(sum(priority_map)) > 0
		tour_times = zeros(length(total_tours), 1);
		for i=1:length(total_tours)
			tour_times(i) = length(total_tours{i});
		end
		[temp, idx] = min(tour_times);
		[ this_energy_cost, this_resource_cost, this_collected_reward, this_wasted_resources, this_tour, priority_map, avoidance_map ] = stochastic_task_GPR( row_distance, vine_distance, priority_map, service_cost_map, beginning, ending, energy_budget, resource_budget, avoidance_map, length(total_tours{idx}));
		total_energy_costs(idx) = total_energy_costs(idx) + this_energy_cost;
		total_resource_costs(idx) = total_resource_costs(idx) + this_resource_cost;
		total_collected_rewards{idx} = [total_collected_rewards{idx}(1:end-1); this_collected_reward(length(total_collected_rewards{idx})+1:end)];
		total_wasted_resources{idx} = [total_wasted_resources{idx}(1:end-1); this_wasted_resources(length(total_wasted_resources{idx})+1:end)];
		total_tours{idx} = [total_tours{idx}(1:end-1); this_tour(length(total_tours{idx})+1:end)];
	end
end
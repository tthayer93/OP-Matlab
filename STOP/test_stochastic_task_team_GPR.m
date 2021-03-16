%test_stochastic_task_team_GPR Run the modified team GPR algorithm on a stochastic task allocation problem
%
%	Version: 1.0
%	Date: 12/17/2020
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This script showcases modified version of the team Greedy Partial Row (GPR) algorithm for the stochastic task allocation problem, on vineyard datasets or randomized datasets, presented in "Task Planning on Stochastic Aisle Graphs for Precision Agriculture" by Kan, Thayer, Carpin, and Karydis.

%% Initialize
clear all;
close all hidden;
%rng(1);
%addpath(genpath('~/Sync/Grad_School/Code'));
%% Parameters
num_rows = 20;
num_vines_per_row = 17;
start_vertex = num_vines_per_row * floor(num_rows/2) + 1;
end_vertex = start_vertex;
num_robots = 2;
energy_budget = 4 * num_vines_per_row;
resource_budget = 2 * num_vines_per_row;
use_vineyard = 1;
create_plots = 1;
if use_vineyard
	%% Use Vineyard data
	desired_moisture = 14;
	vine1_lat =  36.836681;
	vine1_long = -120.214323;
	vineEnd_lat = 36.8432;
	vineEnd_long = -120.209963;
	% filename = 'vineyard_data/07-10-17.csv'; %240x500; 37.491481,-120.539135; 37.484707,-120.530779
	% filename = 'vineyard_data/07-28-17.csv'; %240x500; 37.484536,-120.539205; 37.477375,-120.530805
	dataset = cell(0);
	dataset{end+1} = 'vineyard_data/2018-04-20_ripperdan.csv'; %275x214; 36.836681,-120.214323; 36.8432,-120.209963
	dataset{end+1} = 'vineyard_data/2018-05-07_ripperdan.csv'; %275x214; 36.836681,-120.214323; 36.8432,-120.209963
	dataset{end+1} = 'vineyard_data/2018-06-06_ripperdan.csv'; %275x214; 36.836681,-120.214323; 36.8432,-120.209963
	dataset{end+1} = 'vineyard_data/2018-06-21_ripperdan.csv'; %275x214; 36.836681,-120.214323; 36.8432,-120.209963
	dataset{end+1} = 'vineyard_data/2018-07-05_ripperdan.csv'; %275x214; 36.836681,-120.214323; 36.8432,-120.209963
	dataset{end+1} = 'vineyard_data/2018-07-17_ripperdan.csv'; %275x214; 36.836681,-120.214323; 36.8432,-120.209963
	dataset{end+1} = 'vineyard_data/2018-07-26_ripperdan.csv'; %275x214; 36.836681,-120.214323; 36.8432,-120.209963
	dataset{end+1} = 'vineyard_data/2018-08-06_ripperdan.csv'; %275x214; 36.836681,-120.214323; 36.8432,-120.209963
	dataset{end+1} = 'vineyard_data/2018-08-20_ripperdan.csv'; %275x214; 36.836681,-120.214323; 36.8432,-120.209963
	dataset{end+1} = 'vineyard_data/2018-09-19_ripperdan.csv'; %275x214; 36.836681,-120.214323; 36.8432,-120.209963
	num_trials = length(dataset);
else
	%% Use random data
	priority_multipliers = [1.5, 2];
	num_tasks = 225;
	num_trials = 10;
end
%% Loop over simulated environments
reward_per_visited = [];
wasted_per_visited = [];
visited_vertices = [];
for route_num=1:num_trials
	%% Create priority map and service cost map
	if use_vineyard
		[vertex_table, edge_list, dist_list, rewards, reward_map, row_distance, vine_distance] = build_vineyard_graph(vine1_lat, vine1_long, vineEnd_lat, vineEnd_long, num_rows, num_vines_per_row, dataset{route_num}, 'linear_nearest', 0, desired_moisture, 2, 1);
		start_vertex = (num_vines_per_row+2) * floor(num_rows/2) + 1;
		end_vertex = start_vertex;
		dist_list(:) = 1;
		priority_map = double(reward_map > 0);
		service_cost_map = priority_map .* reward_map;
	else
		priority_map = zeros(num_rows, num_vines_per_row);
		filled = 0;
		while filled < num_tasks
			xy = [randi(num_rows), randi(num_vines_per_row-2)+1];
			if priority_map(xy(1), xy(2)) == 0
				priority_map(xy(1), xy(2)) = priority_multipliers(randi(length(priority_multipliers)));
				filled = filled + 1;
			end
		end
		service_cost_map = double(logical(priority_map)) .* random('Exponential', 2, num_rows, num_vines_per_row);
	end
	row_distance = 1;
	vine_distance = 1;
	%% try team GRP
	try
		[ total_energy_costs, total_resource_costs, total_collected_rewards, total_wasted_resources, total_tours ] = stochastic_task_team_GPR_series( row_distance, vine_distance, priority_map, service_cost_map, start_vertex, end_vertex, energy_budget, resource_budget, num_robots );
% 		avoidance_map = cell(size(priority_map));
% 		total_energy_cost = 0;
% 		total_resource_cost = 0;
% 		total_reward = [];
% 		total_wasted = [];
% 		total_tour = [];
% 		while sum(sum(priority_map)) > 0
% 			[ energy_cost, resource_cost, collected_reward, wasted_resources, tour, priority_map, avoidance_map ] = stochastic_task_GPR( row_distance, vine_distance, priority_map, service_cost_map, start_vertex, end_vertex, energy_budget, resource_budget, avoidance_map, length(total_tour) );
% 			total_energy_cost = total_energy_cost + energy_cost;
% 			total_resource_cost = total_resource_cost + resource_cost;
% 			total_reward = [total_reward(1:end-1); collected_reward(length(total_reward)+1:end)];
% 			total_wasted = [total_wasted(1:end-1); wasted_resources(length(total_wasted)+1:end)];
% 			total_tour = [total_tour(1:end-1); tour(length(total_tour)+1:end)];
% 		end
	catch error
		getReport(error)
		warning("Caught bad route.");
		continue;
	end
	%% Calculate statistics
	all_rewards = [];
	all_wasted = [];
	for i = 1:num_robots
		all_rewards = [all_rewards; total_collected_rewards{i}];
		all_wasted = [all_wasted; total_wasted_resources{i}];
	end
	idxes = (all_rewards > 0) | (all_wasted > 0);
	reward_per_visited = [reward_per_visited, sum(all_rewards(idxes) ./ sum(sum(priority_map .* service_cost_map))) / sum(idxes)];
	wasted_per_visited = [wasted_per_visited, sum(all_wasted(idxes)) / sum(idxes)];
	visited_vertices = [visited_vertices, sum(idxes)];
	%% Plot route results
	if create_plots
		linewidth = 2;
		fontsize = 15;
		y1 = cumsum(all_rewards) ./ sum(sum(priority_map .* service_cost_map));
		y2 = cumsum(all_wasted);
		y1 = y1(idxes);
		y2 = y2(idxes);
		%x = 1:sum(sum(logical(priority_map)));
		x = 1:sum(idxes);
		f1 = figure('visible', 'on');
		hold on;
		set(gcf,'Position',[0, 0, 600, 400]);
		plot(x, y1, 'displayname', 'Greedy Partial Row', 'linewidth', linewidth);
		title('Collected Reward', 'FontSize', fontsize);
		xlabel('Visited vertices');
		ylabel('Percentage of collected reward');
		ax = gca;
		ax.LineWidth = linewidth;
		ax.FontSize = fontsize;
		legend(ax, 'show', 'Location', 'Northwest');
		f2 = figure('visible', 'on');
		hold on;
		set(gcf,'Position',[0, 0, 600, 400]);
		plot(x, y2, 'displayname', 'Greedy Partial Row', 'linewidth', linewidth);
		title('Wasted Resources', 'FontSize', fontsize);
		xlabel('Visited vertices');
		ylabel('Wasted Resources');
		ax = gca;
		ax.LineWidth = linewidth;
		ax.FontSize = fontsize;
		legend(ax, 'show', 'Location', 'Northwest');
	end
	create_plots = 0;
end
%% Print final stats
reward_per_visited
r_v_mean_std = [mean(reward_per_visited), std(reward_per_visited)]
wasted_per_visited
w_v_mean_std = [mean(wasted_per_visited), std(wasted_per_visited)]
visited_vertices
v_v_mean_std = [mean(visited_vertices), std(visited_vertices)]
%RUN_GREEDY_PARTIAL_ROW Run the GPR heuristic on a saved dataset
%
%	Version: 1.0
%	Date: 08/15/2019
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This script showcases the Greedy Row heuristic and Greedy Partial Row (GPR) heuristic presented in https://ieeexplore.ieee.org/abstract/document/8461242

clear all;
%% Load Merced data
% vine1_lat = 37.491481;
% vine1_long = -120.539135;
% vineEnd_lat = 37.484707;
% vineEnd_long = -120.530779;
% num_rows = 240;
% num_vines_per_row = 500;
% filename = 'vineyard_data/07-10-17.csv';
vine1_lat = 37.484536;
vine1_long = -120.539205;
vineEnd_lat = 37.477375;
vineEnd_long = -120.530805;
num_rows = 240;
num_vines_per_row = 500;
filename = 'vineyard_data/07-28-17.csv';
%% Load Ripperdan data
% vine1_lat =  36.836681;
% vine1_long = -120.214323;
% %vineEnd_lat = 36.843812;
% vineEnd_lat = 36.8432;
% vineEnd_long = -120.209963;
% num_rows = 275;
% num_vines_per_row = 214;
% filename = 'vineyard_data/2018-04-20_ripperdan.csv';
% filename = 'vineyard_data/2018-05-07_ripperdan.csv';
% filename = 'vineyard_data/2018-06-06_ripperdan.csv';
% filename = 'vineyard_data/2018-06-21_ripperdan.csv';
% filename = 'vineyard_data/2018-07-05_ripperdan.csv';
% filename = 'vineyard_data/2018-07-17_ripperdan.csv';
% filename = 'vineyard_data/2018-07-26_ripperdan.csv';
% filename = 'vineyard_data/2018-08-06_ripperdan.csv';
% filename = 'vineyard_data/2018-08-20_ripperdan.csv';
% filename = 'vineyard_data/2018-09-19_ripperdan.csv';
%% Setup problem parameters
method = 'linear_nearest';
build_heat_maps = 1;
desired_moisture = 14;
[vertex_table, edge_list, dist_list, rewards, reward_map, row_distance, vine_distance] = build_vineyard_graph(vine1_lat, vine1_long, vineEnd_lat, vineEnd_long, num_rows, num_vines_per_row, filename, method, build_heat_maps, desired_moisture, 1);
beginning = 1;
ending = 1;
starting_budget = num_vines_per_row;
ending_budget = num_vines_per_row*(num_rows+.1*num_rows);
budget_increment = (ending_budget-starting_budget)/15;

%% Normalize edge length to 1
budget = starting_budget:budget_increment:ending_budget;
dist_list(:) = 1;
row_distance = 1;
vine_distance = 1;

%% Solve with heuristics
beginning_row = floor(beginning/num_vines_per_row)+1;
ending_row = floor(ending/num_vines_per_row)+1;
for i=1:length(budget)
    i
	tic;
    [total_cost1(i), total_reward1(i), tour1{i}] = greedy_row_heuristic(row_distance, vine_distance, reward_map, beginning, ending, budget(i));
    time1(i) = toc;
    tic;
    [total_cost2(i), total_reward2(i), tour2{i}] = greedy_partial_row_heuristic(row_distance, vine_distance, reward_map, beginning, ending, budget(i));
    time2(i) = toc;
end
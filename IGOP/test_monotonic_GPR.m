%TEST_MONOTONIC_GPR Test the monotonic GPR algorithm
%
%	Version: 1.0
%	Date: 02/11/2019
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This script runs the monotonic version of GPR, with and without randomization and rolling horizon, to compare against the standard algorithm

clear all;

%% Setup problem parameters
vine1_lat =  36.836681;
vine1_long = -120.214323;
%vineEnd_lat = 36.843812;
vineEnd_lat = 36.8432;
vineEnd_long = -120.209963;
filename = 'vineyard_data/2018-04-20_ripperdan.csv'; %1.0151, 1.0146
filename2 = 'vineyard_data/2018-04-20_ripperdan.mat';
filename3 = '2018-04-20';
% filename = 'Orienteering_Problem/2018-05-07_ripperdan.csv'; %1.01, 1.0074
% filename2 = 'Orienteering_Problem/2018-05-07_ripperdan.mat';
% filename3 = '2018-05-07';
% filename = 'Orienteering_Problem/2018-06-06_ripperdan.csv'; %1.0080, 1.0096
% filename2 = 'Orienteering_Problem/2018-06-06_ripperdan.mat';
% filename3 = '2018-06-06';
% filename = 'Orienteering_Problem/2018-06-21_ripperdan.csv'; %1.0122. 1.0101
% filename2 = 'Orienteering_Problem/2018-06-21_ripperdan.mat';
% filename3 = '2018-06-21';
% filename = 'Orienteering_Problem/2018-07-05_ripperdan.csv'; %1.0082, 1.0076
% filename2 = 'Orienteering_Problem/2018-07-05_ripperdan.mat';
% filename3 = '2018-07-05';
% filename = 'Orienteering_Problem/2018-07-17_ripperdan.csv'; %1.0053, 1.0042
% filename2 = 'Orienteering_Problem/2018-07-17_ripperdan.mat';
% filename3 = '2018-07-17';
% filename = 'Orienteering_Problem/2018-07-26_ripperdan.csv'; %1.0056, 1.0045
% filename2 = 'Orienteering_Problem/2018-07-26_ripperdan.mat';
% filename3 = '2018-07-26';
% filename = 'Orienteering_Problem/2018-08-06_ripperdan.csv'; %1.0064, 1.0066
% filename2 = 'Orienteering_Problem/2018-08-06_ripperdan.mat';
% filename3 = '2018-08-06';
% filename = 'Orienteering_Problem/2018-08-20_ripperdan.csv'; %1.069, 1.0047
% filename2 = 'Orienteering_Problem/2018-08-20_ripperdan.mat';
% filename3 = '2018-08-20';
% filename = 'Orienteering_Problem/2018-09-19_ripperdan.csv'; %1.0030, 1.0030
% filename2 = 'Orienteering_Problem/2018-09-19_ripperdan.mat';
% filename3 = '2018-09-19';

% num_rows = 275;
% num_vines_per_row = 214;
num_rows = 10;
num_vines_per_row = 10;
starting_budget = num_rows;
ending_budget = 1.5*num_rows*num_vines_per_row;
budget_increment = ending_budget/13;
num_agents = 1;
beginning = 1;
ending = 1;
desired_moisture = 13;
min_variance = 0;

best_of = 2;
depth = 8;
num_moves = 1;
deeper = 0;
randomized = 3;
monte_carlo_trials = 32;

%% Load Kriging data
method = 'linear_nearest';
build_heat_maps = 0;
[vertex_table, edge_list, dist_list, rewards, reward_map, row_distance, vine_distance] = build_vineyard_graph(vine1_lat, vine1_long, vineEnd_lat, vineEnd_long, num_rows, num_vines_per_row, filename, method, build_heat_maps, desired_moisture, 1, 0);
coord_list = table2array(vertex_table(:, 2:3));
load(filename2);
moisture_krig = zeros(1, height(vertex_table));
moisture_var = zeros(1, height(vertex_table));
for i=1:height(vertex_table)
    this_x = round(table2array(vertex_table(i, 4)) - xllcorner);
    this_y = round(table2array(vertex_table(i, 5)) - yllcorner);
    moisture_krig(i) = krig_val(this_y, this_x);
    moisture_var(i) = krig_var(this_y, this_x);
end
vertex_table.moisture_krig = moisture_krig';
vertex_table.moisture_var = moisture_var';
variance_map = zeros(num_rows, num_vines_per_row);
v = 0;
for i=1:num_rows
    for j=1:num_vines_per_row
        v = v + 1;
        moisture_map(i, j) = vertex_table.moisture_krig(v);
        this_var = vertex_table.moisture_var(v);
        if this_var >= min_variance
            variance_map(i, j) = this_var;
        end
    end
end
reward_map = abs(desired_moisture - moisture_map);
sampling_map = variance_map;

%% Budget array
budget = (starting_budget:budget_increment:ending_budget);

%% Normalize edge length to 1
dist_list(:) = 1;
row_distance = 1;
vine_distance = 1;

%% Solve with heuristics
tic
for i=1:length(budget)
    [total_cost(i), total_reward(i), tour{i}] = greedy_partial_row_avoidance(vine_distance, row_distance, reward_map, beginning, ending, budget(i));
    [total_cost_1(i), total_reward_1(i), tour_1{i}] = make_vineyard_route_monotonic(vine_distance, row_distance, reward_map, beginning, ending, total_cost(i), total_reward(i), tour{i}, false);
    [total_cost_2(i), total_reward_2(i), tour_2{i}] = insert_partial_rows_to_monotonic_vineyard_route(vine_distance, row_distance, reward_map, budget(i), tour_1{i}, total_cost_1(i), inf);
    [total_cost_3(i), total_reward_3(i), tour_3{i}] = greedy_partial_row_monotonic(vine_distance, row_distance, reward_map, beginning, ending, budget(i));
    [total_cost_4(i), total_reward_4(i), tour_4{i}] = greedy_partial_row_monotonic_horizon_randomized( vine_distance, row_distance, reward_map, beginning, ending, budget(i), depth, best_of, num_moves, deeper, 0);
    [total_cost_5(i), total_reward_5(i), tour_5{i}] = make_vineyard_route_monotonic(vine_distance, row_distance, reward_map, beginning, ending, total_cost_4(i), total_reward_4(i), tour_4{i}, false);
    [total_cost_6(i), total_reward_6(i), tour_6{i}] = insert_partial_rows_to_monotonic_vineyard_route(vine_distance, row_distance, reward_map, budget(i), tour_5{i}, total_cost_5(i), inf);
    [total_cost_7(i), total_reward_7(i), tour_7{i}] = greedy_partial_row_monotonic_horizon_randomized( vine_distance, row_distance, reward_map, beginning, ending, budget(i), depth, best_of, num_moves, deeper, 0, 1);
    total_reward_8(i) = 0;
    total_reward_9(i) = 0;
    total_reward_10(i) = 0;
    for j=1:monte_carlo_trials
        i, j
        [total_cost_j, total_reward_j, tour_j] = greedy_partial_row_monotonic_horizon_randomized( vine_distance, row_distance, reward_map, beginning, ending, budget(i), depth, best_of, num_moves, deeper, randomized);
        if total_reward_j > total_reward_8(i)
            total_cost_8(i) = total_cost_j;
            total_reward_8(i) = total_reward_j;
            tour_8{i} = tour_j;
        end
        [total_cost_k, total_reward_k, tour_k] = make_vineyard_route_monotonic(vine_distance, row_distance, reward_map, beginning, ending, total_cost_j, total_reward_j, tour_j, false);
        if total_reward_k > total_reward_9(i)
            total_cost_9(i) = total_cost_k;
            total_reward_9(i) = total_reward_k;
            tour_9{i} = tour_k;
        end
        [total_cost_l, total_reward_l, tour_l] = insert_partial_rows_to_monotonic_vineyard_route(vine_distance, row_distance, reward_map, budget(i), tour_k, total_cost_k, inf);
        if total_reward_l > total_reward_10(i)
            total_cost_10(i) = total_cost_l;
            total_reward_10(i) = total_reward_l;
            tour_10{i} = tour_l;
        end
    end
end
toc

%% Stats
% percent_diff1 = total_reward_1 ./ max(total_reward_1,  [], 1)
% [~, ~, ranking1] = unique(sum(rank1, 2));
% rank1
% ranking1
% avg_percent_diff1 = sum(percent_diff1, 2)/length(budget)
mean1 = mean(total_reward_1 ./ total_reward, 2)
mean2 = mean(total_reward_2 ./ total_reward, 2)
mean3 = mean(total_reward_3 ./ total_reward, 2)
mean4 = mean(total_reward_4 ./ total_reward, 2)
mean5 = mean(total_reward_5 ./ total_reward, 2)
mean6 = mean(total_reward_6 ./ total_reward, 2)
mean7 = mean(total_reward_7 ./ total_reward, 2)
mean8 = mean(total_reward_8 ./ total_reward, 2)
mean9 = mean(total_reward_9 ./ total_reward, 2)
mean10 = mean(total_reward_10 ./ total_reward, 2)

%% Save Data
wksp_file = sprintf("monotonic_receeding_horizon_heuristics_%s_depth%d_best%d_%dx%db%d-%d-%d.mat", filename3, depth, best_of, num_rows, num_vines_per_row, starting_budget, floor(budget_increment), ending_budget);
save(wksp_file);
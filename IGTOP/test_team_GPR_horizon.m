%TEST_TEAM_GPR_HORIZON Test the team GPR series with rolling horizon algorithm
%
%	Version: 1.0
%	Date: 07/26/2019
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This script runs the team GPR series algorithm using a rolling horizon to compare it against the standard algorithm

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
% filename = 'vineyard_data/2018-05-07_ripperdan.csv'; %1.01, 1.0074
% filename2 = 'vineyard_data/2018-05-07_ripperdan.mat';
% filename3 = '2018-05-07';
% filename = 'vineyard_data/2018-06-06_ripperdan.csv'; %1.0080, 1.0096
% filename2 = 'vineyard_data/2018-06-06_ripperdan.mat';
% filename3 = '2018-06-06';
% filename = 'vineyard_data/2018-06-21_ripperdan.csv'; %1.0122. 1.0101
% filename2 = 'vineyard_data/2018-06-21_ripperdan.mat';
% filename3 = '2018-06-21';
% filename = 'vineyard_data/2018-07-05_ripperdan.csv'; %1.0082, 1.0076
% filename2 = 'vineyard_data/2018-07-05_ripperdan.mat';
% filename3 = '2018-07-05';
% filename = 'vineyard_data/2018-07-17_ripperdan.csv'; %1.0053, 1.0042
% filename2 = 'vineyard_data/2018-07-17_ripperdan.mat';
% filename3 = '2018-07-17';
% filename = 'vineyard_data/2018-07-26_ripperdan.csv'; %1.0056, 1.0045
% filename2 = 'vineyard_data/2018-07-26_ripperdan.mat';
% filename3 = '2018-07-26';
% filename = 'vineyard_data/2018-08-06_ripperdan.csv'; %1.0064, 1.0066
% filename2 = 'vineyard_data/2018-08-06_ripperdan.mat';
% filename3 = '2018-08-06';
% filename = 'vineyard_data/2018-08-20_ripperdan.csv'; %1.069, 1.0047
% filename2 = 'vineyard_data/2018-08-20_ripperdan.mat';
% filename3 = '2018-08-20';
% filename = 'vineyard_data/2018-09-19_ripperdan.csv'; %1.0030, 1.0030
% filename2 = 'vineyard_data/2018-09-19_ripperdan.mat';
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
num_agents = 1:5;
best_of = 5;
depth = 3;

%% Load Kriging data
method = 'linear_nearest';
build_heat_maps = 0;
[vertex_table, edge_list, dist_list, rewards, reward_map, row_distance, vine_distance] = build_vineyard_graph(vine1_lat, vine1_long, vineEnd_lat, vineEnd_long, num_rows, num_vines_per_row, filename, method, build_heat_maps, desired_moisture, 1);
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
avoidance_map = cell(num_rows, num_vines_per_row);
total_cost_2 = zeros(depth, length(budget));
total_reward_2 = zeros(depth, length(budget));
tour_2 = cell(depth, length(budget));
for nAgents = num_agents
    nAgents
    for k=1:depth
        datetime
        k
        for i=1:length(budget)
            if k == 1
                %[total_cost(i), total_reward(i), tour{i}] = greedy_partial_row_avoidance(vine_distance, row_distance, reward_map, beginning, ending, budget(i));
                [total_cost{nAgents}(i, :), total_reward{nAgents}(i), tour{nAgents}{i}] = team_GPR_series(vine_distance, row_distance, reward_map, beginning, ending, budget(i)/nAgents*ones(1, nAgents));
            end
            [ total_cost_2{nAgents}(k, i), total_reward_2{nAgents}(k, i), tour_2{nAgents}{k, i} ] = team_GPR_series_horizon( vine_distance, row_distance, reward_map, beginning, ending, budget(i)/nAgents*ones(1, nAgents), k, best_of, 1, 0, 0);
            total_reward_3{nAgents}(k, i) = 0;
            total_reward_4{nAgents}(k, i) = 0;
            for j=1:10
                [ total_cost_j, total_reward_j, tours_j ] = team_GPR_series_horizon( vine_distance, row_distance, reward_map, beginning, ending, budget(i)/nAgents*ones(1, nAgents), k, best_of, 1, 0, 1);
                if total_reward_j > total_reward_3{nAgents}(k, i)
                    total_cost_3{nAgents}(k, i) = sum(total_cost_j);
                    total_reward_3{nAgents}(k, i) = sum(total_reward_j);
                    tour_3{nAgents}{k, i} = tours_j;
                end
                [ total_cost_j, total_reward_j, tours_j ] = team_GPR_series_horizon( vine_distance, row_distance, reward_map, beginning, ending, budget(i)/nAgents*ones(1, nAgents), k, best_of, k, 0, 1);
                if total_reward_j > total_reward_4{nAgents}(k, i)
                    total_cost_4{nAgents}(k, i) = sum(total_cost_j);
                    total_reward_4{nAgents}(k, i) = sum(total_reward_j);
                    tour_4{nAgents}{k, i} = tours_j;
                end
			end
        end
    end
end
%% Save Data
wksp_file = sprintf("receeding_horizon_heuristics_%s_depth%d_best%d_%dx%db%d-%d-%d.mat", filename3, depth, best_of, num_rows, num_vines_per_row, starting_budget, floor(budget_increment), ending_budget);
save(wksp_file);
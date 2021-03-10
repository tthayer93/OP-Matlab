%% Initialize
clear all;
close all hidden;
rng(2);
addpath(genpath('~/Sync/Grad_School/Code'));
%% Parameters
n_vertex = 10
t_max = 1.25
num_time_steps = 10
%rewards = [[1:n_vertex]', ones(n_vertex, 1)];
%rewards = [[1:n_vertex]', random('exponential', 1, [n_vertex, 1])];
rewards = [[1:n_vertex]', random('uniform', 0, 1, [n_vertex, 1])];
start_vertex = 1;
end_vertex = n_vertex;
alpha = 0.5
cost_dist_type = 'Exponential';
target_failure_rate = 0.1
%% Create vertex set
xy = random('uniform', 0, 1, [n_vertex, 2]);
edge_list = zeros(size(xy, 1)^2, 3);
k = 0;
for i=1:size(xy, 1)
	for j=1:size(xy, 1)
		k = k + 1;
		edge_list(k, :) = [i, j, sqrt((xy(i, 1) - xy(j, 1))^2 + (xy(i, 2) - xy(j, 2))^2)];
	end
end
idxs = (edge_list(:, 1) == edge_list(:, 2));
edge_list(idxs, :) = [];
%alpha = random('uniform', 0, 1, [size(edge_list, 1), 1]);
rewards(start_vertex, 2) = 0;
%% Solve with old method
[op_tour, op_reward, op_cost] = solve_OP(edge_list, rewards(:, 2), t_max, start_vertex, end_vertex);
[vertex_cluster_rewards, vertex_clusters] = OP_route_to_SCOP_expected_costs(op_tour, edge_list, rewards, cost_dist_type, alpha);
[states, state_transition_table, initial_state_distribution, absorbing_states] = clustered_SCOP_to_CMDP_multifail(vertex_cluster_rewards, vertex_clusters, t_max, num_time_steps);
tic
[policy, state_action, rho, old_output] = solve_CMDP(state_transition_table, initial_state_distribution, absorbing_states, target_failure_rate);
toc
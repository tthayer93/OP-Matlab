%% Initialize
clear all;
close all hidden;
%rng(1);
addpath(genpath('~/Sync/Grad_School/Code'));
%% Parameters
n_vertex = 25;
goal_vertexes = n_vertex;
discount_rate = 0.9;
threshold = 0.00001;
max_iterations = 1000;
%% Create table
rewards = zeros(n_vertex, 1);
rewards(goal_vertexes) = 1;
state_transition_table = [];
for i=1:n_vertex
	for j=1:n_vertex
		trans_probs = random('uniform', 0, 1, [n_vertex, 1]);
		trans_probs = trans_probs / sum(trans_probs);
		state_transition_table = [state_transition_table; i*ones(n_vertex, 1), j*ones(n_vertex, 1), [1:n_vertex]', trans_probs, rewards];
	end
end
%% Run value iteration and policy iteration
[values, policy] = value_iteration(state_transition_table, discount_rate, threshold, max_iterations);
[values2, policy2] = policy_iteration(state_transition_table, discount_rate, threshold, max_iterations);
all(policy == policy2)
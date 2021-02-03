%% Initialize
clear all;
close all hidden;
%rng(1);
addpath(genpath('~/Sync/Grad_School/Code'));
%env = rlPredefinedEnv('BasicGridWorld');
env = rlPredefinedEnv('WaterFallGridWorld-Stochastic');
%% Parameters
discount_rate = 0.9;
threshold = 0.001;
max_iterations = 1000;
%% Create table
state_transition_table = [];
for i=1:size(env.Model.T, 1)
	for j=1:size(env.Model.T, 3)
		for k=1:size(env.Model.T, 2)
			state_transition_table = [state_transition_table; i, j, k, env.Model.T(i, k, j), env.Model.R(i, k, j)];
		end
	end
end
%% Run value iteration and policy iteration
[values, policy] = value_iteration(state_transition_table, discount_rate, threshold, max_iterations);
[values2, policy2] = policy_iteration(state_transition_table, discount_rate, threshold, max_iterations);
all(policy == policy2)
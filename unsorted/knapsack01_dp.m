function [opt_val, opt_choice] = knapsack01_dp(values, weights, capacity)
%KNAPSACK01_DP Dynamic programming solver for the 0-1 knapsack problem
%
%	Version: 1.0
%	Date: 14/3/19
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves the 0-1 knapsack problem using dynamic programming.
%	Input:
%		values: a vector containing the rewards for each item
%		weights: a vector containing the weights for each item
%		capacity: the weight capacity of the knapsack
%	Output:
%		opt_val: the total reward of the optimal choice of items
%		opt_choice: a boolean vector determining which items are optimal

    table = zeros(length(weights)+1, capacity+1);
    for i=1:length(weights)
%         for j=1:capacity
%             if weights(i) > j
%                 table(i+1, j+1) = table(i, j+1);
%             else
%                 table(i+1, j+1) = max(table(i, j+1), values(i) + table(i, j-weights(i)+1));
%             end
%         end
        cap = 1:capacity;
        temp = 1:max(length(weights), capacity);
        idxs = temp(weights(i) > cap);
        table(i+1, idxs+1) = table(i, idxs+1);
        idxs = temp(weights(i) <= cap);
        table(i+1, idxs+1) = max(table(i, idxs+1), values(i) + table(i, idxs-weights(i)+1));
    end
    opt_val = table(end,end);
    temp = opt_val;
    i = length(weights); 
    j = capacity;
    opt_choice = zeros(size(values));
    while temp > 0
        while table(i+1,j+1) == temp
            i = i - 1;
        end
        i = i + 1;
        opt_choice(i) = 1;
        j = j - weights(i);
        i = i - 1;
        temp = table(i+1,j+1);
    end

end
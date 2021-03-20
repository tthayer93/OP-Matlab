function [total_cost, total_reward, tour] = make_vineyard_route_monotonic(vine_distance, row_distance, reward_map, beginning, ending, total_cost, total_reward, initial_tour, omit_partials)
%MAKE_VINEYARD_ROUTE_MONOTONIC Eliminate unnecessary movement in a vineyard route
%
%	Version: 1.0
%	Date: 02/05/2019
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function takes in a tour traversing a vineyard and outputs a similar tour that is monotonicly traversing the vineyard.
%	Assumptions:
%		The vineyard is rectangular, such that every row has the same number of vines within it.
%		There are no missing vines within the vineyard (the reward for such vines can be set to 0).
%		Vine rows are equally spaced, and each vine within a row is equally spaced.
%		The agent is allowed to turn around within a row, to access only some of the rewards within it
%		The vine_distance and row_distance are equal to 1 (if this is not the case, using avoidance_map as input or output will not be accurate)
%		Only one robot is allowed within a row at a time
%	Inputs:
%		vine_distance: The distance between each vine in the rows, which is used as the movement cost between them
%		row_distance: The distance between each row of vines, which is used as the movement cost between them
%		reward_map: A matrix of size num_rows*num_vines_per_row containing the reward of each vertex
%		beginning: The vertex at which the tour begins, which must be on the left side (column 1) of the vineyard (reward_map)
%		ending: The vertex at which the tour ends, which must be on the left side (column 1) of the vineyard (reward_map)
%		total_cost: The total movement cost for the computed tour
%		total_reward: The total reward collected for the computed tour
%		initial_tour: A sequence of vertices describing the tour computed by GPR, from beginning vertex to ending vertex
%		omit_partials: A boolean determining whether or not to omit partially traversed rows from the output tour
%	Outputs:
%		total_cost: The total movement cost for the computed tour
%		total_reward: The total reward collected for the computed tour
%		tour: A sequence of vertices describing the newly computed monotonic tour, from beginning vertex to ending vertex
%
% TODO:
% Still needs partials that are not coincident
    
    %% Parameters
    if nargin < 9
        omit_partials = 1;
    end

    %% Extract basic information
    size_row = size(reward_map, 1);
    size_column = size(reward_map, 2);
    vertex_list = 1:size_row*size_column;
    left_side_verticies = vertex_list(mod(vertex_list, size_column) == 1);
    right_side_verticies = vertex_list(mod(vertex_list, size_column) == 0);
    
    %% Figure out which side we are on, left=1 and right=2
    which_side = zeros(length(initial_tour), 1);
    for i=1:length(initial_tour)
        if any(initial_tour(i) == left_side_verticies)
            which_side(i) = 1;
        end
        if any(initial_tour(i) == right_side_verticies)
            which_side(i) = 2;
        end
    end
    
    %% Find partial and full rows
    which_type = zeros(length(initial_tour), 1);
    left_partials = {};
    right_partials = {};
    full_rows = {};
    i = 0;
    while i < length(initial_tour)-1
        i = i + 1;
        if (i > 2) && (which_side(i) ~= 0) && (initial_tour(i+1) == initial_tour(i-1))
            if which_side(i) == 1
                left_partials{end+1} = initial_tour(i);
            else
                right_partials{end+1} = initial_tour(i);
            end
        end
        if (which_side(i) > 0) && (which_side(i+1) == 0)
            this_path = [initial_tour(i)];
            j = 1;
            while which_side(i+j) == 0
                this_path = [this_path; initial_tour(i+j)];
                j = j + 1;
            end
            this_path = [this_path; initial_tour(i+j)];
            if which_side(i) == which_side(i+j)
                for k=0:j
                    which_type(i+k) = 1;
                end
                if which_side(i) == 1
                    left_partials{end+1} = this_path;
                else
                    right_partials{end+1} = this_path;
                end
            else
                for k=0:j
                    which_type(i+k) = 2;
                end
                full_rows{end+1} = this_path;
            end
            i = i + j;
            if i < length(which_side)
                if which_side(i-1) == 0 && which_side(i+1) == 0 && which_side(i) > 0
                    i = i - 1;
                end
            end
        end
    end
    
    %% Sort rows in ascending order (up/down) and remove redundnacy
    left_partials_first = [];
    right_partials_first = [];
    full_rows_first = [];
    i = 0;
    while i < length(left_partials)
        i = i + 1;
        left_partials_first(i) = left_partials{i}(1);
        for j=length(left_partials):-1:(i+1)
            if left_partials{j}(1) == left_partials{i}(1)
                if length(left_partials{j}) > length(left_partials{i})
                    left_partials(i) = [];
                    left_partials_first(i) = [];
                    i = i - 1;
                    break;
                else
                    left_partials(j) = [];
                end
            end
        end
    end
    [temp, idx] = sort(left_partials_first);
    left_partials = {left_partials{idx}};
    left_partials_first = left_partials_first(idx);
    i = 0;
    while i < length(right_partials)
        i = i + 1;
        right_partials_first(i) = right_partials{i}(1);
        for j=length(right_partials):-1:(i+1)
            if right_partials{j}(1) == right_partials{i}(1)
                if length(right_partials{j}) > length(right_partials{i})
                    right_partials(i) = [];
                    right_partials_first(i) = [];
                    i = i - 1;
                    break;
                else
                    right_partials(j) = [];
                end
            end
        end
    end
    [temp, idx] = sort(right_partials_first);
    right_partials = {right_partials{idx}};
    right_partials_first = right_partials_first(idx);
    for i=1:length(full_rows)
        full_rows_first(i) = full_rows{i}(1);
        for j=length(left_partials):-1:1
            if any(full_rows{i} == left_partials{j}(1))
                left_partials(j) = [];
                left_partials_first(j) = [];
            end
        end
        for j=length(right_partials):-1:1
            if any(full_rows{i} == right_partials{j}(1))
                right_partials(j) = [];
                right_partials_first(j) = [];
            end
        end
    end
    [temp, idx] = sort(full_rows_first);
    full_rows = {full_rows{idx}};
    full_rows_first = full_rows_first(idx);
    
    %% Combine monotonicly
    tour = [beginning];
    total_cost = 0;
    which_side = [1];
    while ~isempty(full_rows)
        this_row = ceil(tour(end) / size_column);
        next_row = ceil(full_rows{1}(1) / size_column);
        next_col = mod(full_rows{1}(end) - 1, size_column) + 1;
        direction = 1;
        if tour(end) > full_rows{1}(1)
            direction = -1;
        end
        for i=this_row+direction:direction:next_row
%             if any(tour(end) == left_partials_first)
%                 idx = (tour(end) == left_partials_first);
%                 tour = [tour; left_partials{idx}(2:end)];
%                 total_cost = total_cost + vine_distance*length(left_partials{idx}(2:end));
%                 left_partials(idx) = [];
%                 left_partials_first(idx) = [];
%             end
%             if any(tour(end) == right_partials_first)
%                 idx = (tour(end) == right_partials_first);
%                 tour = [tour; right_partials{idx}(2:end)];
%                 total_cost = total_cost + vine_distance*length(right_partials{idx}(2:end));
%                 right_partials(idx) = [];
%                 right_partials_first(idx) = [];
%             end
            tour = [tour; (i-1)*size_column + which_side(end)];
            total_cost = total_cost + row_distance;
            which_side = [which_side; which_side(end)];
        end
        if next_col == which_side(end)
            full_rows{1} = full_rows{1}(end:-1:1);
        end
        tour = [tour; full_rows{1}(2:end)];
        total_cost = total_cost + vine_distance*length(full_rows{1}(2:end));
        this_side = which_side(end);
        which_side = [which_side; zeros(size(full_rows{1}(2:(end-1))))];
        full_rows(1) = [];
        if this_side == 1
            which_side = [which_side; size_column];
        else
            which_side = [which_side; 1];
        end
    end
    
    %% Insert path to ending vertex
    direction = 1;
    if tour(end) > ending
        direction = -1;
    end
    this_row = ceil(tour(end) / size_column);
    ending_row = ceil(ending / size_column);
    for i=this_row+direction:direction:ending_row
%         if any(tour(end) == left_partials_first)
%             idx = (tour(end) == left_partials_first);
%             tour = [tour; left_partials{idx}(2:end)];
%             total_cost = total_cost + vine_distance*length(left_partials{idx}(2:end));
%             left_partials(idx) = [];
%             left_partials_first(idx) = [];
%         end
        tour = [tour; (i-1)*size_column + 1];
        total_cost = total_cost + row_distance;
        which_side = [which_side; 1];
    end
    
    %% Insert remaining partial rows
    if ~omit_partials
        while ~isempty(left_partials)
            idxs = find(which_side == 1);
            next_row = ceil(left_partials{1}(1) / size_column);
            next_col = mod(left_partials{1}(end) - 1, size_column) + 1;
            distance = abs(next_row - ceil(tour(idxs) / size_column));
            [min_dist, idx] = min(distance);
            prev_row = ceil(tour(idxs(idx)) / size_column);
            subtour = [];
            cost = 0;
            subside = [];
            direction = 1;
            if next_row < prev_row
                direction = -1;
            end
            for i=prev_row+direction:direction:next_row
                subtour = [subtour; (i-1)*size_column + 1];
                cost = cost + row_distance;
                subside = [subside; 1];
            end
            subtour = [subtour; left_partials{1}(2:end)];
            cost = cost + vine_distance*length(left_partials{1}(2:end));
            if ~isempty(subtour)
                if length(subtour) > 1
                    subside = [subside; zeros(size(left_partials{1}(2:end-1))); 1];
                end
            end
            for i=next_row-direction:-direction:prev_row
                subtour = [subtour; (i-1)*size_column + 1];
                cost = cost + row_distance;
                subside = [subside; 1];
            end
            tour = [tour(1:idxs(idx)); subtour; tour(idxs(idx)+1:end)];
            total_cost = total_cost + cost;
            which_side = [which_side(1:idxs(idx)); subside; which_side((idxs(idx)+1):end)];
            left_partials(1) = [];
        end
        while ~isempty(right_partials)
            idxs = find(which_side == size_column);
            next_row = ceil(right_partials{1}(1) / size_column);
            next_col = mod(right_partials{1}(end) - 1, size_column) + 1;
            distance = abs(next_row - ceil(tour(idxs) / size_column));
            [min_dist, idx] = min(distance);
            prev_row = ceil(tour(idxs(idx)) / size_column);
            subtour = [];
            cost = 0;
            subside = [];
            direction = 1;
            if next_row < prev_row
                direction = -1;
            end
            for i=prev_row+direction:direction:next_row
                subtour = [subtour; (i-1)*size_column + size_column];
                cost = cost + row_distance;
                subside = [subside; size_column];
            end
            subtour = [subtour; right_partials{1}(2:end)];
            cost = cost + vine_distance*length(right_partials{1}(2:end));
            if ~isempty(subtour)
                if length(subtour) > 1
                    subside = [subside; zeros(size(right_partials{1}(2:end-1))); size_column];
                end
            end
            for i=next_row-direction:-direction:prev_row
                subtour = [subtour; (i-1)*size_column + size_column];
                cost = cost + row_distance;
                subside = [subside; size_column];
            end
            tour = [tour(1:idxs(idx)); subtour; tour(idxs(idx)+1:end)];
            total_cost = total_cost + cost;
            which_side = [which_side(1:idxs(idx)); subside; which_side((idxs(idx)+1):end)];
            right_partials(1) = [];
        end
    end
    
    %% Calculate reward
    transverse_reward_map = reward_map';
    reward_vector = transverse_reward_map(:);
    total_reward = sum(reward_vector(unique(tour)));
    
end


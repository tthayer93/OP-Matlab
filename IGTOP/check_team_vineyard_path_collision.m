function [ maybe ] = check_team_vineyard_path_collision( reward_map, paths )
%CHECK_TEAM_VINEYARD_PATH_COLLISION
%
%	Version: 1.0
%	Date: 02/17/2020
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function checks for potential robot collisions for the Irrigation Graph Team Orienteering Problem (IGTOP) discussed in https://ieeexplore.ieee.org/abstract/document/9062300
%	Assumptions:
%		The vineyard is rectangular, such that every row has the same number of vines within it.
%		There are no missing vines within the vineyard (the reward for such vines can be set to 0).
%		Vine rows are equally spaced, and each vine within a row is equally spaced.
%		The agent is allowed to turn around within a row, to access only some of the rewards within it
%		The vine_distance and row_distance are equal to 1 (if this is not the case, the results of this check function will be innacurate)
%		Only one robot is allowed within a row at a time
%	Inputs:
%		reward_map: A matrix of size num_rows*num_vines_per_row containing the reward of each vertex
%		paths: A cell array containing the computed paths for each robot, with individual cells containing a single robot's path
%	Outputs:
%		maybe: A boolean determining if there is a collision or not (collision -> maybe = true)

    rows = size(reward_map, 1);
    cols = size(reward_map, 2);
    boundries = [];
    k = 1;
    for i=1:rows
        for j=1:cols
            if j==1 || j==cols
                boundries = [boundries; k];
            end
            k = k + 1;
        end
    end
    maybe = false;
    for i=1:(length(paths)-1)
        for j=(i+1):length(paths)
            for k=2:min(length(paths{i}), length(paths{j}))
                if paths{i}(k-1) > 0 && paths{j}(k-1) > 0 && paths{i}(k) > 0 && paths{j}(k) > 0
                    first_i = paths{i}(k-1);
                    first_j = paths{j}(k-1);
                    second_i = paths{i}(k);
                    second_j = paths{j}(k);
                    %if ~(any(first_i == boundries) && any(first_j == boundries) && any(second_i == boundries) && any(second_j == boundries))
                    if ~((any(first_i == boundries) || any(first_j == boundries)) && (any(second_i == boundries) || any(second_j == boundries)))
                        if ((first_i == second_j) && (second_i == first_j)) || ((first_i == first_j) && (second_i == second_j))
                            maybe = true;
                            %return
                        end
					end
                end
            end
        end
    end

end


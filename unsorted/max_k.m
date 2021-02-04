function [ vals, indexes ] = max_k( arr, k, stable )
%MAX_K Returns the k largest elements of an array
%
%	Version: 1.0
%	Date: unknown
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function quickly (when k << n) finds the largest k elements of an array. O(n*k)
%	Inputs:
%		arr: an array
%		k: the number of values to return
%		stable: a boolean determining whether to return the values in the order they appear in arr (1) or to sort them in descending order (0). default = 1
%	Output:
%		vals: an array of the k largest values in arr
%		indexes: an array of the indices where the values in vals appear in arr 

    if nargin == 2
        stable = 1;
    end
    arr_size = size(arr);
    if k <= length(arr) && k > 0
        vals = arr(1:k);
        indexes = 1:k;
        [min_val, min_idx] = min(vals, [], 'includenan');
        for i=(k+1):length(arr)
            if arr(i) > min_val || isnan(min_val)
                vals(min_idx) = [];
                indexes(min_idx) = [];
                if arr_size(1) < arr_size(2)
                    vals = [vals, arr(i)];
                else
                    vals = [vals; arr(i)];
                end
                indexes = [indexes, i];
                [min_val, min_idx] = min(vals, [], 'includenan');
            end
        end
    elseif k <= 0
        vals = [];
        indexes = [];
    else
        vals = arr;
        indexes = 1:length(arr);
    end
    if ~stable
        [h, idxes] = sort(vals, 'descend');
        vals = vals(idxes);
        indexes = indexes(idxes);
    end

end


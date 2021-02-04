function unique_a = unique_int(a, if_logical) 
%%UNIQUE_INT Fast way to check for unique values in integer array
%
%	Version: 0.0
%	Date: unknown
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function finds the unique values in an integer array quickly.
%	Inputs:
%		a: an array that must contain only positive integers
%		if_logical: a binary variable that determines the format of the output. default = 0
%	Output:
%		unique_a: the unique values in a, with two possible formats
%			if_logical == 0: a vector containing the each unique value of integer in a
%			if_logical == 1: a logical vector of length max(a) with 1 where integer values exist in a and 0 otherwise

    if nargin < 2
        if_logical = 0;
	end
    temp = zeros(1, max(a));
    temp(a) = temp(a) + 1;
    if if_logical
        unique_a = temp;
    else
        unique_a = find(temp);
    end

end
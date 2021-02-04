function C = fastintersect(A,B)
%FASTINTERSECT Quickly find common values in two integer arrays
%
%	Version: 1.0
%	Date: unknown
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	A fast way to find the intersection (common values) of two integer arrays.
%	Inputs:
%		A: an integer array
%		B: an integer array
%	Outputs:
%		C: an integer array of values that occur in both A and B

    if ~isempty(A)&&~isempty(B)
       temp = zeros(1, max(max(A),max(B)) ) ;
       temp(A) = 1;
       C = B(logical(temp(B)));
    else
        C = [];
    end
end
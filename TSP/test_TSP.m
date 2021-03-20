% TEST_TSP Test the ILP solver for an instance of the TSP
%
%	Version: 1.0
%	Date: 07/18/2020
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This script runs the Integer Linear Program (ILP) solver for the Traveling Salesman Problem (TSP) on either a USA map or a random set of vertices.

%% Initialize
clear all;
close all hidden;
%rng(1);
addpath(genpath('~/Sync/Grad_School/Code'));
warning('off','all');
%% Parameters
n_vertex = 25
symmetric = 0;
use_usa = 0;
%% Load data for continental USA
if use_usa
	load('usborder.mat','x','y','xx','yy');
	rng(3,'twister') % makes a plot with stops in Maine & Florida, and is reproducible
	nStops = n_vertex; % you can use any number, but the problem size scales as N^2
	stopsLon = zeros(nStops,1); % allocate x-coordinates of nStops
	stopsLat = stopsLon; % allocate y-coordinates
	n = 1;
	while (n <= nStops)
		xp = rand*1.5;
		yp = rand;
		if inpolygon(xp,yp,x,y) % test if inside the border
			stopsLon(n) = xp;
			stopsLat(n) = yp;
			n = n+1;
		end
	end
	plot(x,y,'Color','red'); % draw the outside border
	hold on
	% Add the stops to the map
	plot(stopsLon,stopsLat,'*b')
	hold off
	xy = [stopsLon, stopsLat];
else
	xy = random('uniform', 0, 1, [n_vertex, 2]);
end
%% Create edge list
edge_list = zeros(size(xy, 1)^2, 3);
k = 0;
for i=1:size(xy, 1)
	if symmetric
		for j=i:size(xy, 1)
			k = k + 1;
			edge_list(k, :) = [i, j, sqrt((xy(i, 1) - xy(j, 1))^2 + (xy(i, 2) - xy(j, 2))^2)];
		end
	else
		for j=1:size(xy, 1)
			k = k + 1;
			edge_list(k, :) = [i, j, sqrt((xy(i, 1) - xy(j, 1))^2 + (xy(i, 2) - xy(j, 2))^2)];
		end
	end
end
idxs = (edge_list(:, 1) == edge_list(:, 2));
edge_list(idxs, :) = [];
%% Solve with intlinprog
tic
[tour, total_cost] = solve_TSP(edge_list, symmetric);
toc
%% Plot
hold on;
scatter(xy(:, 1), xy(:, 2), '*');
for i=1:length(tour)-1
	quiver(xy(tour(i), 1),xy(tour(i), 2),xy(tour(i+1), 1)-xy(tour(i), 1),xy(tour(i+1), 2)-xy(tour(i), 2),1,'color',[0,0,1]);
end
quiver(xy(tour(end), 1),xy(tour(end), 2),xy(tour(1), 1)-xy(tour(end), 1),xy(tour(1), 2)-xy(tour(end), 2),1,'color',[0,0,1]);

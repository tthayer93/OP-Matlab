function [ vertex_table, edge_list, dist_list, rewards, reward_map, row_distance, vine_distance ] = build_vineyard_graph( vine1_lat, vine1_long, vineEnd_lat, vineEnd_long, num_rows, num_vines_per_row, filename, method, build_heat_maps, desired_moisture )
%BUILD_VINEYARD_GRAPH Builds a graph using vineyard data
%
%	Version: 1.0
%	Date: 09/03/2017
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function takes soil moisture data from a vineyard and interpolates it to create a graph with vertex rewards.
%	Assumptions:
%		Reward values are calculated as R = abs(desired_moisture - interpolated_moisture)
%		The vineyard is rectangular, such that every row has the same number of vines within it.
%		There are no missing vines within the vineyard (the reward for such vines can be set to 0 after this function returns).
%		Vine rows are equally spaced, and each vine within a row is equally spaced.
%	Inputs:
%		vine1_lat: Latitude of the upper left vine
%		vine1_long: Longitude of the upper left vine
%		vineEnd_lat: Latitude of the lower right vine
%		vine_End_long: Longitude of the lower right vine
%		num_rows: Number of vine rows in the vineyard
%		num_vines_per_row: Number of vines per row in the vineyard
%		filename: The filename and path for the csv file in which the moisture data is contained, which must have these columns
%			'TimeStamp','Record','Zone','Latitude','Longitude','Moisture','Periodus','Attenuation','Permittivity','ProbeModel'
%		method: The interpolation method for calculating moisture values of vines not directly measured, which may be either 'nearest_nearest' or 'linear_nearest'
%		build_heat_maps: A boolean defining whether to build a visual heat map of the vineyard
%		desired_moisture: A constant for the desired moisture level of the entire vineyard
%	Outputs:
%		vertex_table: A table of vertices in the graph containing the following data
%			Vertex | x | y | utm_N | utm_E | utm_Zone | Latitude | Longitude | Moisture | Rewards
%		edge_list: A list of valid edges within the graph, with the following format
%			start_vertex | end_vertex
%		dist_list: A list of distances or cost for each edge, with each row corresponding to the edge in the same row of edge_list
%		rewards: A vector containing the reward for each vertex in the graph
%		reward_map: A matrix of size num_rows*num_vines_per_row containing the reward of each vertex
%		row_distance: The distance between each row of vines
%		vine_distance: The distance between each vine in the rows

    %% Define spacing of vines if needed
    [vine1_utm_n, vine1_utm_e, vine1_utm_zone] = deg2utm(vine1_lat, vine1_long);
    [vineEnd_utm_n, vineEnd_utm_e, vineEnd_utm_zone] = deg2utm(vineEnd_lat, vineEnd_long);
    if vine1_utm_zone ~= vineEnd_utm_zone
        'UTM Zones do not match. Ending script...'
        return;
    end
    row_dist = (vineEnd_utm_e - vine1_utm_e) / num_rows;
    vine_dist = (vineEnd_utm_n - vine1_utm_n) / num_vines_per_row;

    %% Build xy map
    nVertex = num_rows * num_vines_per_row;
    vertex = [1:nVertex]';
    x = zeros(size(vertex));
    y = zeros(size(vertex));
    for i=1:num_rows
        for j=1:num_vines_per_row
            n = num_vines_per_row*(i-1)+j;
            x((i-1)*num_vines_per_row+j) = (j-1)*vine_dist; %x
            y((i-1)*num_vines_per_row+j) = (i-1)*row_dist; %y
        end
    end
    [vine1_utm_n, vine1_utm_e, vine1_utm_zone] = deg2utm(vine1_lat, vine1_long);
    utm_N = x + vine1_utm_n;
    utm_E = y + vine1_utm_e;
    utm_Zone = repmat(vine1_utm_zone, size(utm_E));
    [latitude, longitude] = utm2deg(utm_N, utm_E, utm_Zone);
    vertex_table = table(vertex, x, y, utm_N, utm_E, utm_Zone, latitude, longitude,'VariableNames', {'Vertex', 'x', 'y', 'utm_N', 'utm_E', 'utm_Zone', 'Latitude', 'Longitude'});

    %% Build edge_list and dist_list
    edge_list = [];
    dist_list = [];
    for i=1:num_rows
        if i ~= num_rows
            edge_list(end+1,:) = [num_vines_per_row*(i-1)+1, num_vines_per_row*(i)+1];
            dist_list(end+1,:) = row_dist;
            edge_list(end+1,:) = [num_vines_per_row*(i)+1, num_vines_per_row*(i-1)+1];
            dist_list(end+1,:) = row_dist;
        end
        for j=2:num_vines_per_row
            edge_list(end+1,:) = [num_vines_per_row*(i-1)+j-1, num_vines_per_row*(i-1)+j];
            dist_list(end+1,:) = vine_dist;
            edge_list(end+1,:) = [num_vines_per_row*(i-1)+j, num_vines_per_row*(i-1)+j-1];
            dist_list(end+1,:) = vine_dist;
        end
        if i ~= num_rows
            edge_list(end+1,:) = [num_vines_per_row*(i), num_vines_per_row*(i+1)];
            dist_list(end+1,:) = row_dist;
            edge_list(end+1,:) = [num_vines_per_row*(i+1), num_vines_per_row*(i)];
            dist_list(end+1,:) = row_dist;
        end
    end
    dist_list = abs(dist_list);

    %% Build moisture map estimate
    fileID = fopen(filename,'r');
    try
        fileID = fopen(filename,'r');
        format_spec = '%{MM/dd/yyyy hh:mm:ss a}D%f%C%f%f%f%f%f%f%C%[^\n\r]';
        dataArray = textscan(fileID, format_spec, 'Delimiter', ',', 'TextType', 'string', 'HeaderLines' , 1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    catch
        fileID = fopen(filename,'r');
        format_spec = '%{MM/dd/yyyy HH:mm}D%f%q%f%f%f%f%f%f%q%[^\n\r]';
        dataArray = textscan(fileID, format_spec, 'Delimiter', ',', 'TextType', 'string', 'HeaderLines' , 1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    end
    fclose(fileID);
    moisture_data = table(dataArray{1:end-1}, 'VariableNames', {'TimeStamp','Record','Zone','Latitude','Longitude','Moisture','Periodus','Attenuation','Permittivity','ProbeModel'});
    moisture_funct = scatteredInterpolant(table2array(moisture_data(:,[4, 5])), table2array(moisture_data(:, 6)));
    if strcmp(method, 'nearest_nearest')
        moisture_funct.Method = 'nearest'; %nearest linear natural
        moisture_funct.ExtrapolationMethod = 'nearest';
    elseif strcmp(method, 'linear_nearest')
        moisture_funct.Method = 'linear';
        moisture_funct.ExtrapolationMethod = 'nearest';
    else
        "Only 'nearest_nearest' and 'linear_nearest' are valid methods. Defaulting to linear_nearest."
        moisture_funct.Method = 'linear';
        moisture_funct.ExtrapolationMethod = 'nearest';
    end
    vertex_table.Moisture = moisture_funct(vertex_table.Latitude, vertex_table.Longitude);

    %% Build rewards
    vertex_table.Rewards = abs(desired_moisture - vertex_table.Moisture);
    rewards = vertex_table.Rewards';

    %% Build heat maps
    moisture_map = zeros(num_rows, num_vines_per_row);
    reward_map = zeros(num_rows, num_vines_per_row);
    v = 0;
    for i=1:num_rows
        for j=1:num_vines_per_row
            v = v + 1;
            moisture_map(i, j) = vertex_table.Moisture(v);
            reward_map(i, j) = vertex_table.Rewards(v);
        end
    end
    if build_heat_maps == 1
        hold on;
        figure(1);
        colormap('cool');
        imagesc(moisture_map);
        colorbar;
        figure(2);
        colormap('jet');
        imagesc(reward_map);
        colorbar;
    end

    %% Compile rewards and costs for each row
    row_distance = abs(row_dist);
    vine_distance = abs(vine_dist);

end
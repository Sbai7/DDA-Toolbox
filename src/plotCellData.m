function plotCellData(g, data, varargin)
%plotCellData Visualize scalar or index-based data on a Cartesian grid.
%
%   plotCellData(g, data) visualizes cell-centered data on a 2D or 3D
%   Cartesian grid. Each cell of the grid is colored based on values in 
%   'data', which can either represent scalar values or class indices.
%
%   plotCellData(g, data, 'DataType', dataType) specifies the type of data:
%       - 'Scalar': (default) treats 'data' as scalar values and uses a 
%         continuous colormap based on the data range.
%       - 'Index': treats 'data' as a set of discrete class indices.
%         Each unique index in 'data' is assigned a distinct color.
%
%   plotCellData(g, data, 'DataType', dataType, Name, Value) allows for
%   additional Name-Value pair arguments to control appearance. These
%   options are passed to 'plotGrid' and subsequently to MATLAB's 'patch'
%   function for rendering.
%
%   INPUT PARAMETERS:
%   g       - Struct representing a Cartesian grid, such as one created
%             by the 'cartGrid' function. The grid structure 'g' must
%             include:
%                - g.ne: the number of cells (elements) in the grid.
%   data      - Array containing either scalar data or class indices.
%               For 'DataType' set to 'Scalar', data values must be numeric 
%               scalars. For 'DataType' set to 'Index', data values must 
%               be finite integers representing class indices.
%   varargin  - Optional Name-Value pairs, including:
%               - 'DataType': 'Scalar' or 'Index' (default: 'Scalar').
%
%   Name-Value Pair Arguments:
%   These are optional arguments that control properties such as color, 
%   transparency, or lighting, which will be applied to the grid plot.
%
%   Example:
%       % Scalar data example
%       g = cartGrid([10, 10]);        % 2D grid of 10x10 cells
%       data = rand(g.ne, 1);          % Random scalar data per cell
%       plotCellData(g, data);         % Visualize scalar data on the grid
%
%       % Index data example
%       data = randi([1, 4], g.ne, 1); % Random indices representing classes
%       plotCellData(g, data, 'DataType', 'Index'); % Visualize class indices
%
%   NOTE: The function displays a colorbar that represents the range of 
%   scalar data values across the grid.
%
%   Author: M.A. Sbai, Ph.D.
%
% Copyright (C) 2024 Mohammed Adil Sbai
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.


    % Separate 'DataType' parameter from other Name-Value pairs
    dataType = 'Scalar'; % Default value for DataType
    remainingArgs = {};  % Collect other Name-Value pairs

    % Parse varargin manually to allow flexible Name-Value pairs
    for i = 1:2:length(varargin)
        if strcmp(varargin{i}, 'DataType')
            dataType = varargin{i+1};
        else
            remainingArgs = [remainingArgs, varargin{i}, varargin{i+1}];
        end
    end

    % Validate data size consistency with the grid structure
    assert(g.ne == size(data, 1), ...
        'The number of rows in data must match the number of grid cells.');

    % Retrieve the current colormap and setup variables
    cmap = colormap();             % Get current colormap
    nc   = size(cmap, 1);          % Number of colors in colormap

    % Check DataType and process data accordingly
    if strcmp(dataType, 'Index')
        % 'Index' mode: data is a set of discrete class indices
        unique_classes = unique(data(:));      % Find unique class indices
        num_classes = numel(unique_classes);   % Count distinct classes

        % Select colors corresponding to the number of classes
        color_idx = round(linspace(1, nc, num_classes));
        colors = cmap(color_idx, :);           % Map each class to a distinct color

        % Initialize color_data with appropriate size
        color_data = zeros(g.ne, 3);

        % Assign colors to each class
        for i = 1:num_classes
            class_mask = data == unique_classes(i);    % Logical mask for class
            color_data(class_mask, :) = repmat(colors(i, :), sum(class_mask), 1);
        end

    else
        % 'Scalar' mode: data is continuous and normalized to the colormap
        data_min = min(data); 
        data_max = max(data);
        rng = [data_min, data_max];            % Data range for normalization

        % Handle edge case where all data values are the same
        if rng(2) <= rng(1)
            color_data = repmat(cmap(ceil(nc / 2), :), g.ne, 1);
        else
            % Normalize data to fit within colormap
            ix = ceil(nc * (data(:) - rng(1)) / diff(rng));
            ix = max(1, min(ix, nc));          % Ensure indices are in bounds
            color_data = cmap(ix, :);          % Map data to colors
        end
    end

    % Plot the grid with cell colors based on the computed color data
    plotGrid(g, 'FaceColor', 'flat', 'CData', reshape(color_data, g.ne, 1, 3), remainingArgs{:});

    % Set color axis and add colorbar based on data type
    if strcmp(dataType, 'Index')
        % Index mode: discrete color bar labels
        clim([1, num_classes]);                % Set limits for class indices
        colormap(colors);                      % Use selected colors for classes
        cb = colorbar;
        cb.Ticks = linspace(1, num_classes, num_classes);
        cb.TickLabels = arrayfun(@num2str, unique_classes, 'UniformOutput', false);
    else
        % Scalar mode: continuous color bar
        set(gca, 'CLim', [data_min, data_max]);
        colorbar;
    end
end
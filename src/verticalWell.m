function W = verticalWell(W, G, i, j, k_layers, rate, varargin)
% Add a vertical well to an existing well structure for finite difference grids
%
% function W = verticalWell(G, W, i, j, k_layers, rate, varargin)
%
% Add a well to the set of wells in the structure W, positioned at grid indices 
% (i, j, k_layers) and with a total injection or pumping rate specified in 'rate'.
% The 'Type' field must be either 'injection' or 'pumping'.
% 
% Arguments:
% - G: Grid structure (should contain the fields Nx, Ny, Nz representing grid dimensions).
% - W: Existing well structure (can be empty if this is the first well).
% - i, j: Grid indices corresponding to the (i, j) position of the well in the grid.
% - k_layers: A vector specifying the vertical grid layers the well spans.
% - rate: Positive scalar value for the total flow rate (injection or pumping),
%         a function handle for a transient rate, or a time series of time-rate pairs.
% - varargin: Optional arguments for additional well properties (e.g., 'radius', 'name').
%
% Examples:
%   W = verticalWell([], G, 10, 20, [1, 2, 3], @(t) 300 * exp(-0.1 * t), 'type', 'injection', 'name', 'Well1', 'radius', 0.1);
%   W = verticalWell(W, G, 5, 15, [2], [0 200; 1 250; 2 300; 4 150], 'type', 'pumping', 'name', 'Well2', 'radius', 0.15);
%
% Author: M.A. Sbai, Ph.D.
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


% Default well radius and type
radius = 0.1; % Default radius if not provided
wellName = ''; % Default well name
wellType = ''; % Initialize well type

% Parse optional arguments
for iVar = 1:2:numel(varargin)
    switch lower(varargin{iVar})
        case 'type'
            wellType = varargin{iVar+1};
            % Ensure Type is either 'injection' or 'pumping'
            if ~ismember(wellType, {'injection', 'pumping'})
                error('Well Type must be either ''injection'' or ''pumping''.');
            end
        case 'radius'
            radius = varargin{iVar+1};
        case 'name'
            wellName = varargin{iVar+1};
        otherwise
            error(['Unknown property: ', varargin{iVar}]);
    end
end

% Ensure the 'rate' value is appropriate
if isscalar(rate) && rate <= 0
    error('Well ''rate'' must be a positive value.');
elseif ~isscalar(rate) && ~isa(rate, 'function_handle') && ~ismatrix(rate)
    error('Well ''rate'' must be a positive scalar, a function handle, or a time series of time-rate pairs.');
elseif ismatrix(rate) && ~isscalar(rate)
    % Check if rate is a time-rate pair matrix
    if size(rate, 2) ~= 2 || any(rate(:, 1) < 0) || any(rate(:, 2) <= 0)
        error('Time-rate pairs must be in the form [time rate], with positive rates.');
    end
end

% Check that the grid indices are within the grid dimensions
if i > G.Nx || j > G.Ny || any(k_layers > G.Nz)
    error('Well indices exceed the grid dimensions.');
end

% Determine the grid cells corresponding to the well's position (i, j, k_layers)
cells = sub2ind([G.Nx, G.Ny, G.Nz], repmat(i, size(k_layers)), repmat(j, size(k_layers)), k_layers);

% Append the new well information to the well structure
new_well.i = i;
new_well.j = j;
new_well.k_layers = k_layers;
new_well.radius = radius;
new_well.cells = cells;
new_well.type = wellType; % 'injection' or 'pumping'
new_well.name = wellName; % Name of the well
new_well.rate = rate;

% Append the new well to the existing well structure
W = [W, new_well];

end

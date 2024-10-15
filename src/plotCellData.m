function plotCellData(g, data, varargin)
%plotCellData: Plots cell-centered scalar data on a given 2D/3D cartesian 
%              grid.
%
% INPUTS:
% g                 - Cartesian grid structure as created by calling 
%                     cartGrid routine   
% data              - array of scalar data to visualize on external facets
%                     of the grid 
% varargin (optional)
%                   - optional property-value pair arguments submitted to 
%                     PlotGrid routine and subsequently to MATLAB patch 
%                     command. Type 'help patch' in the command window to 
%                     see the list of supported options
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



assert(g.ne==size(data,1), ...
    'Number of rows in data must equal the number of cells in the grid!');

% save data min/max
data_min = min(data); 
data_max = max(data); 

% transform scalar data to RGB values 
cmap = colormap();
rng  = [min(data), max(data)];
nc   = size(cmap, 1);

if ~(rng(2) > rng(1))
    data = cmap(ceil(nc / 2), :);
else
	ix   = ceil(nc * (data(:) - rng(1)) ./ diff(rng));
	data = cmap(max(1, min(ix, nc)), :);
end

plotGrid(g,'FaceColor','flat','CData',reshape(data,g.ne,1,3),varargin{:});

set(gca,'CLim',[data_min data_max]), colorbar

end
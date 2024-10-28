function [wellPair, wellPairIdx, wellPairVolume] = wellPairCatchments(Grid, ...
    injWellIdx, pumpWellIdx)
%wellPairCatchments Calculates aquifer catchments owned by well pairs.
%   
% INPUTS:
%   injWellIdx  - N-by-1 grid cell array of injection well indices as 
%                 returned from injectionCatchments function
%   pumpWellIdx - N-by-1 grid cell array of pumping well indices as
%                 returned from pumpingCatchments function
%
% OUTPUTS:
%   wellPair    - Np-by-2 array of injection-pumping well pairs that own 
%                 a catchment in the aquifer. Np is the total number of
%                 well pairs. The first column is injection wells with 'I'
%                 appended, and the second column is pumping wells with 'P'
%                 appended.
%   wellPairIdx - N-by-1 grid cell array of well pair indices. That is, each
%                 cell identifies to which injection-pumping well pair it 
%                 belongs.
%   wellPairVolume  - Np-by-1 array of aquifer pore volumes of all identified
%                 injection-pumping well pairs in the aquifer system.
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

% Get the number of grid cells
N = Grid.N;

% Check that inputs vectors have the same size
assert(N == size(injWellIdx,1) & N == size(pumpWellIdx,1) & ...
    size(injWellIdx,2) == 1 & size(pumpWellIdx,2) == 1, ... 
    'Input arrays size are not correct.');

% Find unique injection-pumping well pairs, starting indices of each well
% pair, and an array of grid cell well pair indices. 
[wellPairNum,~,wellPairIdx] = unique([injWellIdx,pumpWellIdx], ... 
    'rows','stable'); 

% Initialize wellPair as a string array
wellPair = strings(size(wellPairNum, 1), 2);

% Append 'I' to injection well indices (first column) and 'P' to pumping
% well indices (second column)
for i = 1:size(wellPairNum, 1)
    wellPair(i, 1) = ['I' num2str(wellPairNum(i, 1))];  % 'I' + injection well index
    wellPair(i, 2) = ['P' num2str(wellPairNum(i, 2))];  % 'P' + pumping well index
end

% Calculate aquifer pore volume in each injection-pumping well pair.
wellPairVolume = zeros(size(wellPairNum,1),1);
for i = 1 : size(wellPairNum, 1)
   %num_cells = nnz(wellPairIdx == i);
   cells_porosity = Grid.properties.porosity(wellPairIdx == i);
   wellPairVolume(i) = sum(cells_porosity) * Grid.V;
end

end

function at = travelTime(G, V, W, forward)
% Calculates gridded arrival times of a tracer in forward or backward
% directions.
%
% INPUTS:
%   G    - Grid structure used for discretization 
%   V    - Structure of mid-cells edges velocities (x, y, z)
%   W    - Wells structure containing wells positions, rates, names, etc.
% forward- Logical flag indicating flow direction (true for forward, false for backward)
%
% OUTPUTS:
% at     - Array of cell-centered travel times for a chemical tracer
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


% Validate inputs
assert(islogical(forward), 'Input "forward" must be a logical value.');
N = G.N; 

% Get the wells count
if forward
   wellIndices = find(strcmp({W.type}, 'injection'));
else 
   wellIndices = find(strcmp({W.type}, 'pumping'));
end

% Get the number of involved wells 
numWells = length(wellIndices);

% If no injection or pumping wells are present, return empty array
if numWells == 0
    at = [];
    return;
end

% Get well flow rates (Qw) array
Qw = getFlowRateVector(N, W, 0);

% Calculate arrival times from injection or production wells
at = steadyTravelTime(forward, G, V, Qw);

end
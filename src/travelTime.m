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
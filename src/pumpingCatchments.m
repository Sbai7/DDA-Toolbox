function [partition, volume, nd_indicator] = pumpingCatchments(G, V, W)
%pumpingCatchments: Determines the capture zones of all pumping wells,
% computes the capture volumes, and evaluates the numerical dispersion.
%
% INPUTS:
%   G    - Structure containing grid discretization details
%   V    - Specific discharge field (structure with x, y, z components)
%   W    - Wells structure containing wells positions, rates, names, etc.
%  
% OUTPUTS:
%   partition    - Grid array giving the index of the pumping well to which 
%                  each grid cell belongs. 
%   volume       - Array of capture volumes for each pumping well
%   nd_indicator - Array of numerical dispersion indicators
%
% Author: M.A. Sbai, Ph.D.

% Number of grid cells
N = G.N; 

% Reverse velocity field
V.x = -V.x; 
V.y = -V.y; 
V.z = -V.z; 

% Identify and count pumping wells
pumpingIndices = find(strcmp({W.type}, 'pumping'));
numPumpingWells = length(pumpingIndices);

% Extract the pumping wells
pumpingWells = W(pumpingIndices);

% exit if there are any wells!
if numPumpingWells == 0 
    partition = zeros(N,1);
    volume = zeros(numPumpingWells,1);
    nd_indicator = zeros(numPumpingWells,1);
    return; 
end

% Initialize well flow rates (Qw) array
Qw = getFlowRateVector(N, W, 0);

% Calculate capture zone for each pumping well
C_wells = zeros(N, numPumpingWells);
for w = 1:numPumpingWells
    C_ext = zeros(N, 1);                           % initialize external conc.
    C_ext(pumpingWells(w).cells) = 1;              % fixed unit concentration 
    C = steadyConcentration(G, V, C_ext, -Qw);     % solve for steady-state transport
    C_wells(:, w) = C;                             % assign this well contribution
end

% Majority vote over tracer concentrations
[~, partition] = max(C_wells, [], 2);

% Set cells with no capture to 0
partition(all(C_wells == 0, 2)) = 0;

% Calculate pore volumes for each pumping well capture region
volume = zeros(numPumpingWells, 1);
por = G.properties.porosity(:);
for w = 1:numPumpingWells
    p = find(partition == w);
    volume(w) = sum(G.V .* por(p));
end

% Evaluate numerical dispersion indicator for each capture zone
nd_indicator = zeros(numPumpingWells, 1);
for w = 1:numPumpingWells
    idx = find(C_wells(:, w));
    nd_indicator(w) = 100 * norm(C_wells(idx, w) - ones(numel(idx), 1)) / numel(idx);
end

end
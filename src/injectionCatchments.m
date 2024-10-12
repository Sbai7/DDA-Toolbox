function [partition, volume, nd_indicator] = injectionCatchments(G, V, W)
%injectionCatchments: Calculates the catchment zones for injection wells.
%
% INPUTS:
%   G    - Grid structure used for discretization 
%   V    - Structure of mid-cells edges velocities (x, y, z)
%   W    - Wells structure containing wells positions, rates, names, etc.
%
% OUTPUTS: 
%   partition   - Gridded array giving the injection well index to which
%                 each grid cell belongs.
%   volume      - Array of captured volumes by each injection well
%   nd_indicator- Numerical dispersion indicator
%
% Author: M.A. Sbai, Ph.D.
%

% Number of grid cells
N = G.N; 

% Identify and count injection wells
injectionIndices = find(strcmp({W.type}, 'injection'));
numInjectionWells = length(injectionIndices);

% Extract the injection wells
injectionWells = W(injectionIndices);

% exit if there are any wells!
if numInjectionWells == 0
   partition = zeros(N,1);
   volume  = zeros(numInjectionWells,1);
   nd_indicator = zeros(numInjectionWells,1);
   return;
end

%--- Initialize well flow rates (Qw) array
Qw = getFlowRateVector(N, W, 0);

% Calculate drainage zone for each injection well
C_wells = zeros(N, numInjectionWells);
for w = 1:numInjectionWells
    C_ext = zeros(N, 1);
    C_ext(injectionWells(w).cells) = 1;
    C = steadyConcentration(G, V, C_ext, Qw);    
    C_wells(:, w) = C;
end

% Majority vote over tracer concentrations to determine partition
[~, partition] = max(C_wells, [], 2); 

% set cells with no drainage to 0
partition(all(C_wells == 0, 2)) = 0;

% Calculate pore volumes of each injection well swept region 
volume = zeros(numInjectionWells, 1);
por = G.properties.porosity(:);
for w = 1:numInjectionWells
    p = find(partition == w);
    volume(w) = sum(G.V .* por(p));
end

% Evaluate numerical dispersion indicator for each drainage zone 
nd_indicator = zeros(numInjectionWells, 1);
for w = 1:numInjectionWells
    idx = find(C_wells(:, w));
    nd_indicator(w) = 100 * norm(C_wells(idx, w) - ones(numel(idx), 1)) / numel(idx);
end

end
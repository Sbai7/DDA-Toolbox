function [h, V] = head(Grid, W, h, time, dt)
% head - Solves the 3D hydraulic or groundwater head equation in  
% steady-state and transient modes.
%
% Syntax: [h, V] = head(Grid, W, h, time, dt)
%
% Inputs:
%   Grid - Structure containing grid discretization details
%   W    - Wells structure containing wells positions, rates, names, etc.
%   h (optional)  - Head from the previous time step (transient mode)
%   time (optional) - Current simulation time
%   dt (optional) - Time interval for transient solver (transient mode)
%
% Outputs:
%   h - Array of cell-centered heads
%   V - Structure containing Darcian velocities along x, y, and z 
%       directions at mid-cells
%
% Author: M.A. Sbai, Ph.D.

%--- Initialize well flow rates (Qw) array
if nargin == 2
   Qw = getFlowRateVector(Grid.N, W, 0);
else
   Qw = getFlowRateVector(Grid.N, W, time);
end

%--- Compute hydraulic head and specific discharge fields
if nargin > 4
    % Call finite-difference solver in transient mode
    [h, V] = tpfa(Grid, Qw, h, dt);
else
    % Call finite-difference solver in steady-state mode
    [h, V] = tpfa(Grid, Qw);
end

end

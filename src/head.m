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

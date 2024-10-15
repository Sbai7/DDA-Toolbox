function tau = steadyTravelTime(forward, Grid, V, Qw)
% Calculates cell-centered tracer arrival times in either forward or backward
% directions using a stationary, fully implicit first-order upwind method.
%
% INPUTS:
% forward           - Logical flag indicating the direction of flow (true for forward, false for backward)
% Grid              - Grid used for discretization 
% V                 - Structure containing mid-cell edge velocities in x, y, and z directions
% Qw                - Array of injection/production flow rates 
%
% OUTPUT:
% tau               - Array of cell-centered travel times
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

% Compute the upwind advection matrix
if forward 
    A = upwindAdvectionMatrix(Grid, V, Qw);
else
    % Reverse velocity field and use -Qw for backward direction
    V.x = -V.x; V.y = -V.y; V.z = -V.z;
    A = upwindAdvectionMatrix(Grid, V, -Qw);
    % Restore the original velocity field
    V.x = -V.x; V.y = -V.y; V.z = -V.z;
end

% initialize linear solver if not given as input argument
% ls       = LinearSolver();
% ls.maxit = 1000;
% ls.tol   = 1e-12;
% ls.ilutype = 'nofill'; 
% tau = ls.solve(-A,Grid.V(:).*Grid.por(:),'Solver','bicgstab');
% fprintf('%s\n',ls.convergenceMessage());
% 

% Compute travel times using MATLAB's backslash operator for solving linear systems
tau = -A \ (Grid.V(:) .* Grid.properties.porosity(:));

end
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
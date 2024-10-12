function C = steadyConcentration(Grid, V, externalConcentration, qw, ls)
%steadyConcentration: Solves the 3D convective part of the mass transport
% equation using a fully implicit first-order upwind scheme.
%
% INPUTS:
%   Grid              - Grid structure used for discretization
%   V                 - Structure of mid-cells edges velocities (x, y, z)
%   externalConcentration - Array of external concentrations (injected or inflowing)
%   qw                - Injection/Production flow rates
%   ls (optional)     - Linear solver and preconditioner types for solving the system
%
% OUTPUTS:
%   C                 - Array of cell-centered steady-state concentrations
%
% Author: M.A. Sbai, Ph.D.
%

% Number of unknowns
N = Grid.Nx * Grid.Ny * Grid.Nz;

% Compute Upwind stencil sparse system matrix
A = upwindAdvectionMatrix(Grid, V, qw);

% Initialize linear solver if not provided
if nargin < 5
    ls = LinearSolver();
    ls.maxit = 1000;
    ls.tol = 1e-6;
    ls.ilutype = 'ilutp'; % Preconditioner type
end
        
% Concentration at injection cells
injectionConcentration = Grid.V(:) .* max(qw, 0) .* externalConcentration;

% Solve linear system -A*C = injectionConcentration using iterative solver
C = ls.solve(-A, injectionConcentration, 'Solver', 'bicgstab');
fprintf('%s\n', ls.convergenceMessage());
    
%     C = callAMGCL(speye(N)-B,C0+Ci, 'solver', 'bicgstab', ...
%                     'preconditioner', 'amg',...
%                     'coarsening', 'ruge_stuben',...                  
%                     'relaxation', 'gauss_seidel', ...
%                     'maxIterations', 10000, ...
%                     'verbose', true);
end
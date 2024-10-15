function C = steadyConcentration(Grid, V, externalConcentration, Qw, ls)
%steadyConcentration: Solves the 3D convective part of the mass transport
% equation using a fully implicit first-order upwind scheme.
%
% INPUTS:
%   Grid              - Grid structure used for discretization
%   V                 - Structure of mid-cells edges velocities (x, y, z)
%   externalConcentration - Array of external concentrations (injected or inflowing)
%   Qw                - Injection/Production flow rates
%   ls (optional)     - Linear solver and preconditioner types for solving the system
%
% OUTPUTS:
%   C                 - Array of cell-centered steady-state concentrations
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


% Number of unknowns
N = Grid.Nx * Grid.Ny * Grid.Nz;

% Compute Upwind stencil sparse system matrix
A = upwindAdvectionMatrix(Grid, V, Qw);

% Initialize linear solver if not provided
if nargin < 5
    ls = LinearSolver();
    ls.maxit = 1000;
    ls.tol = 1e-6;
    ls.ilutype = 'ilutp'; % Preconditioner type
end
        
% Concentration at injection cells
injectionConcentration = Grid.V(:) .* max(Qw, 0) .* externalConcentration;

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
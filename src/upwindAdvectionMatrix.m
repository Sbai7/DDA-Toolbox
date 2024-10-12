function A = upwindAdvectionMatrix(Grid, V, q)
%UpwindAdvectionMatrix: Generates the pure upwind advection matrix for a 
% mass/heat/particle transport process based on the spatial distribution of 
% the velocity field.
%
% INPUTS:
%   Grid - Grid structure used for discretization
%   V    - Structure of mid-cells edges velocities (x, y, z)
%   q    - Injection/Production flow rates
%
% OUTPUT:
%   A    - Sparse matrix of the advection process
% 
% Author: M.A. Sbai, Ph.D.
% 

% Number of grid cells
Nx = Grid.Nx; 
Ny = Grid.Ny; 
Nz = Grid.Nz; 
N = Grid.N;

% Extract production flow rates (i.e., negative part)
productionFlowRates = min(q, 0);

% Separate velocity components into positive and negative parts
xNegative = min(V.x, 0);
yNegative = min(V.y, 0);
zNegative = min(V.z, 0);
xPositive = max(V.x, 0);
yPositive = max(V.y, 0);
zPositive = max(V.z, 0);

% Reshape arrays into column vectors
xNegative = reshape(xNegative(1:Nx,:,:), N, 1);
yNegative = reshape(yNegative(:,1:Ny,:), N, 1);
zNegative = reshape(zNegative(:,:,1:Nz), N, 1);
xPositive = reshape(xPositive(2:Nx+1,:,:), N, 1);
yPositive = reshape(yPositive(:,2:Ny+1,:), N, 1);
zPositive = reshape(zPositive(:,:,2:Nz+1), N, 1);

% Upwind advection sparse matrix stencil
% Diagonal vectors and indices
diagonalVectors = [zPositive, yPositive, xPositive, ...
                   productionFlowRates + xNegative - xPositive + yNegative - yPositive + zNegative - zPositive, ...
                   -xNegative, -yNegative, -zNegative];
diagonalIndices = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];

% Construct sparse matrix
A = spdiags(diagonalVectors, diagonalIndices, N, N);

end
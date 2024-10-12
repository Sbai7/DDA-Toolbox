function [h, V] = tpfa(Grid, Qw, h, dt)
% tpfa - Solves the 3D Groundwater equation using Two-Point Flux  
%        Approximation under steady-state or transient conditions.
%
% Syntax: [h, V] = tpfa(Grid, Qw, h, dt)
%
% Inputs:
%   Grid - Structure containing grid discretization details
%   Qw   - A N-by-1 array of cellwise injection/production flow rates
%   h (optional)  - Head from the previous time step (transient mode)
%   dt (optional) - Time interval for transient solver (transient mode)
%
% Outputs:
%   h - Array of cell-centered heads
%   V - Structure containing Darcian velocities along x, y, and z 
%       directions at mid-cells
%
% Author: M.A. Sbai, Ph.D.

Nx = Grid.Nx; 
Ny = Grid.Ny; 
Nz = Grid.Nz; 
N = Grid.N;

hx = Grid.hx; 
hy = Grid.hy; 
hz = Grid.hz;

% Compute transmissibilities using harmonic averaging
L = Grid.properties.hydraulic_conductivity.^(-1);  % inverse of the hydraulic conductivity tensor
tx = 2 * hy * hz / hx; 
ty = 2 * hx * hz / hy; 
tz = 2 * hx * hy / hz; 

TX = zeros(Nx+1, Ny, Nz); 
TY = zeros(Nx, Ny+1, Nz); 
TZ = zeros(Nx, Ny, Nz+1);

% Compute transmissibility values
TX(2:Nx, :, :) = tx ./ (L(1, 1:Nx-1, :, :) + L(1, 2:Nx, :, :));
TY(:, 2:Ny, :) = ty ./ (L(2, :, 1:Ny-1, :) + L(2, :, 2:Ny, :));
TZ(:, :, 2:Nz) = tz ./ (L(3, :, :, 1:Nz-1) + L(3, :, :, 2:Nz));

% Assemble the TPFA discretization matrix
x1 = reshape(TX(1:Nx, :, :), N, 1); 
x2 = reshape(TX(2:Nx+1, :, :), N, 1);
y1 = reshape(TY(:, 1:Ny, :), N, 1); 
y2 = reshape(TY(:, 2:Ny+1, :), N, 1);
z1 = reshape(TZ(:, :, 1:Nz), N, 1); 
z2 = reshape(TZ(:, :, 2:Nz+1), N, 1);

DiagVecs = [-z2, -y2, -x2, x1 + x2 + y1 + y2 + z1 + z2, -x1, -y1, -z1];
DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
A = spdiags(DiagVecs, DiagIndx, N, N);

% Handle transient conditions if necessary
if nargin == 4
    stor = Grid.specific_storage / dt;
    tterm = reshape(stor, N, 1);
    A = A + spdiags(tterm, 0, N, N);
    h = reshape(h, N, 1);
    Qw = Qw + (tterm .* h);
end

%--- Solve resulting sparse linear system of equations using AMGCL
[h, err, nIter] = callAMGCL(A, Qw, 'solver', 'cg', ...
                    'preconditioner', 'amg',...
                    'coarsening', 'aggregation',...                  
                    'relaxation', 'ilu0', ...
                    'maxIterations', 10000, ...
                    'tolerance',     1e-06);
fprintf('Number of PCG iterations = %d\n', nIter);
fprintf('Residual norm error      = %G\n', err);

% Solve the sparse linear system
% ls = LinearSolver();
% ls.maxit = 10000; 
% ls.tol = 1e-12;
% P = ls.solve(A, q, 'Solver', 'pcg');
% fprintf('%s\n', ls.convergenceMessage());
%ls.plot();
%saveas(gcf, ['bench1_conv_',ls.icholtype,'_',num2str(Nx),'x',num2str(Ny),'x',num2str(Nz)], 'tiff');

h = reshape(h, Nx, Ny, Nz);

% Calculate mid-cells edges fluxes using Darcy's law
V.x = zeros(Nx+1, Ny, Nz); 
V.y = zeros(Nx, Ny+1, Nz); 
V.z = zeros(Nx, Ny, Nz+1);

V.x(2:Nx, :, :) = (h(1:Nx-1, :, :) - h(2:Nx, :, :)) .* TX(2:Nx, :, :);
V.y(:, 2:Ny, :) = (h(:, 1:Ny-1, :) - h(:, 2:Ny, :)) .* TY(:, 2:Ny, :);
V.z(:, :, 2:Nz) = (h(:, :, 1:Nz-1) - h(:, :, 2:Nz)) .* TZ(:, :, 2:Nz);

end
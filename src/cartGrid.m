function G = cartGrid(x, y, z)
% Create a uniform/non-uniform Cartesian grid from rectilinear arrays along
% each dimension.
%
% G = cartGrid(x)
% G = cartGrid(x, y)
% G = cartGrid(x, y, z)
%
% Examples:
% 1) Create a 2D Cartesian uniform grid over the domain 0 < x < 1000 and 0
% < y < 2000 with a 10 and 20m spacing, respectively:
% >> G = cartGrid(0:10:1000, 0:20:2000); 
% >> plotGrid(G, 'FaceColor', 'none', 'EdgeColor', 'r');
%
% 2) Create a 3D uniform grid in X and Y, but non-uniform in Z over the
% domain 0 < x < 1000, 0 < y < 2000, and -50 < z < 25:
% >> G = cartGrid(0:10:1000, 0:20:2000, [25 20 15 10 5 0 -5 -15 -30 -50]);
% >> plotGrid(G, 'FaceColor', 'c', 'FaceAlpha', 0.7);
%
% Author: M.A. Sbai, Ph.D.

% Validate input dimensions and ensure non-decreasing grid coordinates
if nargin == 1
    y = 1; % Default 1D case
    z = 1;
elseif nargin == 2
    z = 1; % Default 2D case
end

% Error checking for non-decreasing coordinates
if any(diff(x) <= 0), error('x-coordinates must be strictly increasing.'); end
if nargin > 1 && any(diff(y) <= 0), error('y-coordinates must be strictly increasing.'); end
if nargin > 2 && any(diff(z) <= 0), error('z-coordinates must be strictly increasing.'); end

% Get grid dimensions
G.ni = length(x);
G.nj = length(y);
G.nk = length(z);
G.nn = G.ni * G.nj * G.nk;

% Total number of elements
if nargin == 3
    G.ne = (G.ni-1) * (G.nj-1) * (G.nk-1); % 3D case
elseif nargin == 2
    G.ne = (G.ni-1) * (G.nj-1); % 2D case
else
    G.ne = G.ni-1; % 1D case
end

% Dimension and orientation
G.dimension = nargin;
G.orientation = sprintf('%dD', G.dimension);

% Create node coordinates for the grid
if nargin == 3
    [X, Y, Z] = ndgrid(x, y, z);
    G.coord = [X(:), Y(:), Z(:)];
elseif nargin == 2
    [X, Y] = ndgrid(x, y);
    G.coord = [X(:), Y(:)];
else
    G.coord = x(:);
end

% Preallocate element connectivity based on the dimension
if nargin == 3
    G.elem_nodes = zeros(G.ne, 8); % 8 nodes per element in 3D
elseif nargin == 2
    G.elem_nodes = zeros(G.ne, 4); % 4 nodes per element in 2D
else
    G.elem_nodes = zeros(G.ne, 2); % 2 nodes per element in 1D
end

% Calculate element-to-node connectivity in a vectorized manner
if nargin == 3
    % 3D element connectivity
    [I, J, K] = ndgrid(1:G.ni-1, 1:G.nj-1, 1:G.nk-1);
    idx = @(i, j, k) i + (j-1) * G.ni + (k-1) * G.ni * G.nj;
    ie = (1:G.ne)';
    
    G.elem_nodes(ie,1) = idx(I(:), J(:), K(:));
    G.elem_nodes(ie,2) = idx(I(:)+1, J(:), K(:));
    G.elem_nodes(ie,3) = idx(I(:)+1, J(:)+1, K(:));
    G.elem_nodes(ie,4) = idx(I(:), J(:)+1, K(:));
    G.elem_nodes(ie,5) = idx(I(:), J(:), K(:)+1);
    G.elem_nodes(ie,6) = idx(I(:)+1, J(:), K(:)+1);
    G.elem_nodes(ie,7) = idx(I(:)+1, J(:)+1, K(:)+1);
    G.elem_nodes(ie,8) = idx(I(:), J(:)+1, K(:)+1);
    
elseif nargin == 2
    % 2D element connectivity
    [I, J] = ndgrid(1:G.ni-1, 1:G.nj-1);
    ie = (1:G.ne)';
    
    G.elem_nodes(ie,1) = I(:) + (J(:)-1) * G.ni;
    G.elem_nodes(ie,2) = G.elem_nodes(ie,1) + 1;
    G.elem_nodes(ie,3) = G.elem_nodes(ie,1) + G.ni;
    G.elem_nodes(ie,4) = G.elem_nodes(ie,3) + 1;
    
else
    % 1D element connectivity
    G.elem_nodes(:,1) = (1:G.ne)';
    G.elem_nodes(:,2) = G.elem_nodes(:,1) + 1;
end

% Optional: Identify boundary nodes (useful for boundary conditions)
G.boundary_nodes = identifyBoundaryNodes(G);

% Create the Grid structure and assign properties
G.Nx = G.ni - 1;
G.Ny = G.nj - 1;
G.Nz = G.nk - 1;
G.N  = G.ne;

% Calculate grid spacing (cell size in each direction)
G.hx = (x(end) - x(1)) / G.Nx;
G.hy = (y(end) - y(1)) / G.Ny;
G.hz = (z(end) - z(1)) / G.Nz;

% Calculate the volume of each grid cell
G.V = G.hx * G.hy * G.hz;

end


function bnodes = identifyBoundaryNodes(G)
% Identify boundary nodes for a given grid G
% This function returns a logical array that marks which nodes are on the
% boundary of the domain.

bnodes = false(G.nn, 1);

if G.dimension == 1
    bnodes(1) = true;
    bnodes(end) = true;

elseif G.dimension == 2
    % Mark the boundary nodes for a 2D grid
    bnodes(1:G.ni) = true; % Bottom
    bnodes(end-G.ni+1:end) = true; % Top
    bnodes(1:G.ni:G.nn) = true; % Left
    bnodes(G.ni:G.ni:G.nn) = true; % Right

elseif G.dimension == 3
    % Mark the boundary nodes for a 3D grid
    ni = G.ni; nj = G.nj; nk = G.nk;
    
    % X-Y faces at Z=constant
    bnodes(1:ni*nj) = true; % Z = min(z)
    bnodes(end-ni*nj+1:end) = true; % Z = max(z)
    
    % Y-Z faces at X=constant
    bnodes(1:ni:ni*nj*nk) = true; % X = min(x)
    bnodes(ni:ni:ni*nj*nk) = true; % X = max(x)
    
    % X-Z faces at Y=constant
    bnodes(reshape(1:ni*nj:nk*ni*nj, [], 1)) = true; % Y = min(y)
    bnodes(reshape(ni:ni:ni*nj:nk*ni*nj, [], 1)) = true; % Y = max(y)
end


end

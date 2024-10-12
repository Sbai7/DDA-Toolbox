function hs = plotStreamline(V, h, s, forward, varargin)
% Plots streamlines traced by the velocity structure returned by the 
% Pressure solver. 
%
% INPUTS:
% V                 - Structure of mid-cells edges velocities along x,y,and
%                     z directions
% h                 - Uniform spacings of the grid along X,Y & Z directions
% s                 - 2D or 3D array of starting streamlines X,Y, & Z 
%                     positions. The array size is 2-by-ns or 3-by-ns where
%                     ns is the number of streamlines to trace. The leading
%                     dimension indicates the space dimension (2 or 3). 
% forward           - Flag for forward (true) or backward (false) tracing
%                     of streamlines
%
% OUTPUTS:
% hs                - Graphics handle to plotted streamlines by this
%                     function
%
%
% Author: M.A. Sbai, Ph.D.
%             

dim = size(s,1);

assert(dim==2 || dim==3);
assert(length(h)==dim);
assert(isa(forward,'logical'));

Nx = size(V.x,1)-1;
Ny = size(V.y,2)-1;

if dim==2

    Ux = [zeros(1,Ny); V.x(2:Nx,:); zeros(1,Ny)]; 
    Ux = 0.5*(Ux(1:end-1,:)+Ux(2:end,:));
    Uy = [zeros(Nx,1), V.y(:,2:Ny), zeros(Nx,1)];
    Uy = 0.5*(Uy(:,1:end-1)+Uy(:,2:end));

    hx = h(1); hy = h(2);
    [X,Y] = meshgrid((1:Nx)*hx-0.5*hx,(1:Ny)*hy-0.5*hy);
    if ~forward, Ux=-Ux; Uy=-Uy; end    
    hs    = streamline(X,Y,Ux',Uy',s(1,:),s(2,:), varargin{:}) ;

else
    
    error('PlotStreamline does not support 3D streamlines yet!');
    
end

end


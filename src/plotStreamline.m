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


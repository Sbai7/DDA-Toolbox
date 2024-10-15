function Qw = getFlowRateVector(N, W, time)
%getFlowRateVector Returns a N-by-1 cell-wise array of all pumping or
%injection flow rates.
%
% Syntax: Qw = getFlowRateVector(N, W, time)
%
% Inputs:
%   N    - Number of grid cells
%   W    - Wells structure containing wells positions, rates, names, etc.
%   time - Current simulation time
%
% Outputs:
%   Qw   - A N-by-1 array of cell-wise injection/production flow rates
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
Qw = zeros(N, 1);

%--- Distribute the flow rate across well cells
for w = 1 : numel(W)
   
   % Check for constant flow rate case
   if isscalar(W(w).rate) && isnumeric(W(w).rate)
       cell_flow_rate = W(w).rate / numel(W(w).cells);
   elseif isa(W(w).rate, 'function_handle')
       % Evaluate the anonymous function for time-varying rate
       cell_flow_rate = W(w).rate(time) / numel(W(w).cells);
   elseif ismatrix(W(w).rate) && size(W(w).rate, 2) == 2
       % Interpolate from the time series
       % W(w).rate is expected to be a [n x 2] matrix where first column is time, second is rate
       cell_flow_rate = interp1(W(w).rate(:, 1), W(w).rate(:, 2), time, 'makima', 'extrap') / numel(W(w).cells);
   else
       error('Unrecognized rate format for well %s', W(w).Name);
   end

   % Make sure that pumping wells have negative flow rates
   if strcmpi(W(w).type, 'pumping')
       cell_flow_rate = -cell_flow_rate;
   elseif ~strcmpi(W(w).type, 'injection')
       error('Invalid well type for well %s. Type must be "injection" or "pumping".', W(w).Name);
   end
   
   % Assign the flow rate to all well cells
   Qw(W(w).cells) = Qw(W(w).cells) + cell_flow_rate;
end

end
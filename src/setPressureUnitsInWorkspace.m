function setPressureUnitsInWorkspace()
% setPressureUnitsInWorkspace: Sets pressure variables (Pascal, Bar, Atmosphere, Torr)
% directly in the base workspace.
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

    % Define the pressure units in Pascals
    pascal = 1;                        % 1 Pa
    bar = 1e5;                         % 1 bar = 100,000 Pa
    atmosphere = 101325;               % 1 atm = 101,325 Pa
    torr = 133.322;                    % 1 Torr = 133.322 Pa
    psi = 6894.76;                     % 1 psi = 6894.76 Pa

    % Assign variables to the base workspace
    assignin('base', 'pascal', pascal);
    assignin('base', 'bar', bar);
    assignin('base', 'atmosphere', atmosphere);
    assignin('base', 'torr', torr);
    assignin('base', 'psi', psi);
end

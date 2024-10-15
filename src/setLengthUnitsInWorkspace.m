function setLengthUnitsInWorkspace()
% setLengthUnitsInWorkspace: Sets length variables (meter, kilometer, centimeter, millimeter)
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

    % Define the length units in meters
    meter = 1;                   % 1 meter
    kilometer = 1000 * meter;    % 1 kilometer = 1000 meters
    centimeter = 0.01 * meter;    % 1 centimeter = 0.01 meters
    millimeter = 0.001 * meter;   % 1 millimeter = 0.001 meters
    micrometer = 1e-6 * meter;    % 1 micrometer = 1e-6 meters
    nanometer = 1e-9 * meter;      % 1 nanometer = 1e-9 meters
    mile = 1609.34 * meter;        % 1 mile = 1609.34 meters
    yard = 0.9144 * meter;         % 1 yard = 0.9144 meters
    foot = 0.3048 * meter;         % 1 foot = 0.3048 meters
    inch = 0.0254 * meter;         % 1 inch = 0.0254 meters

    % Assign variables to the base workspace
    assignin('base', 'meter', meter);
    assignin('base', 'kilometer', kilometer);
    assignin('base', 'centimeter', centimeter);
    assignin('base', 'millimeter', millimeter);
    assignin('base', 'micrometer', micrometer);
    assignin('base', 'nanometer', nanometer);
    assignin('base', 'mile', mile);
    assignin('base', 'yard', yard);
    assignin('base', 'foot', foot);
    assignin('base', 'inch', inch);
end

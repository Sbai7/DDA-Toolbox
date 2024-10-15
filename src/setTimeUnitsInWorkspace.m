function setTimeUnitsInWorkspace(year_value)
% setTimeUnitsInWorkspace: Sets time variables (second, minute, hour, day, week, month, year)
% directly in the base workspace.
%
% Inputs:
%   year_value - (Optional) Specify the year to check if it's a leap year.
%                If no year is provided, defaults to 365 days for a non-leap year.
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

    
    if nargin < 1
        year_value = 365;  % Default to 365 days if no input is given.
    end

    % Define the time units
    second = 1;
    minute = 60 * second;
    hour   = 60 * minute;
    day    = 24 * hour;
    week   = 7 * day;       % Define a week
    month  = 30 * day;      % Approximate average month

    % Leap year detection and year calculation
    if isLeapYear(year_value)
        year = 366 * day;   % Leap year has 366 days
    else
        year = 365 * day;   % Non-leap year has 365 days
    end

    % Assign variables to the base workspace
    assignin('base', 'second', second);
    assignin('base', 'minute', minute);
    assignin('base', 'hour', hour);
    assignin('base', 'day', day);
    assignin('base', 'week', week);
    assignin('base', 'month', month);
    assignin('base', 'year', year);
end


function leap = isLeapYear(year_value)
    % isLeapYear: Checks if a given year is a leap year.
    %
    % Inputs:
    %   year_value - The year to check for leap year.
    %
    % Output:
    %   leap - Boolean value, true if it's a leap year, false otherwise.

    if mod(year_value, 4) == 0
        if mod(year_value, 100) == 0
            if mod(year_value, 400) == 0
                leap = true;  % Divisible by 400
            else
                leap = false; % Divisible by 100 but not by 400
            end
        else
            leap = true;      % Divisible by 4 but not by 100
        end
    else
        leap = false;         % Not divisible by 4
    end
end

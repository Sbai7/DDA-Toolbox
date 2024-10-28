function intersectionPoint = findIntersectionX(contourData, contourIndex, x0)
    % findIntersectionX Find the intersection of a contour with a vertical line at x = x0.
    %
    % Syntax:
    %   intersectionPoint = findIntersectionX(contourData, contourIndex, x0)
    %
    % Description:
    %   This function finds the intersection point (x0, y) of a specified contour line
    %   with a vertical line at x = x0, using linear interpolation between contour
    %   vertices. The function returns the coordinates of the intersection point or
    %   an empty array if no intersection is found.
    %
    % Input Arguments:
    %   contourData   - (cell array) A cell array of structures where each structure
    %                   contains the fields:
    %                     - x: x-coordinates of the contour line vertices.
    %                     - y: y-coordinates of the contour line vertices.
    %   contourIndex  - (integer scalar) Index of the contour line in `contourData`
    %                   to search for an intersection.
    %   x0            - (numeric scalar) The x-coordinate of the vertical line with
    %                   which to find the intersection.
    %
    % Output:
    %   intersectionPoint - (1-by-2 numeric array) The coordinates [x0, y] of the
    %                       intersection point, where y is interpolated. Returns an
    %                       empty array if no intersection is found.
    %
    % Example:
    %   % Given a contourData structure
    %   intersection = findIntersectionX(contourData, 1, 3.5);
    %
    % Notes:
    %   - If the specified contour index is invalid, an error is thrown.
    %   - If no intersection is found, the function returns an empty array and
    %     displays a warning.
    %
    % See also: struct, warning
    
    % Check if the contour index is valid
    if contourIndex < 1 || contourIndex > length(contourData)
        error('Invalid contour index.');
    end
    
    % Retrieve the x and y coordinates of the specified contour
    xCoords = contourData{contourIndex}.x;
    yCoords = contourData{contourIndex}.y;
    
    % Initialize the intersection point
    intersectionPoint = []; 
    
    % Loop through the vertices and check for intersections
    for i = 1:length(xCoords)-1
        % Check if x0 is between the x values of the current segment
        if (xCoords(i) <= x0 && xCoords(i+1) >= x0) || (xCoords(i) >= x0 && xCoords(i+1) <= x0)
            % Linear interpolation to find the corresponding y value
            % Calculate the slope of the line segment
            slope = (yCoords(i+1) - yCoords(i)) / (xCoords(i+1) - xCoords(i));
            % Interpolate to find the y value at x0
            yIntersect = yCoords(i) + slope * (x0 - xCoords(i));
            intersectionPoint = [x0, yIntersect]; % Store the intersection point
            return; % Exit the function after finding the first intersection
        end
    end
    
    % If no intersection is found, return an empty array
    if isempty(intersectionPoint)
        warning('No intersection found for the specified contour at x = %.2f.', x0);
    end
end

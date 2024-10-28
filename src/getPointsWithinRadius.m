function points = getPointsWithinRadius(xy, r, numPoints)
    % getPointsWithinRadius Generate points within a radius around given coordinates.
    %
    % Syntax:
    %   points = getPointsWithinRadius(xy, r, numPoints)
    %
    % Description:
    %   This function generates a set of points around each coordinate in
    %   the input array `xy` within a specified radius `r`. It returns an 
    %   array containing `numPoints` evenly distributed points around each 
    %   center point within the circle defined by `r`.
    %
    % Input Arguments:
    %   xy        - (n-by-2 numeric array) Each row of `xy` represents a
    %               point's (x, y) coordinates.
    %   r         - (numeric scalar) Radius of the circle within which the
    %               points are generated around each input coordinate.
    %   numPoints - (integer scalar) Number of points to generate around 
    %               each center point, distributed evenly on the circle.
    %
    % Output:
    %   points - (2-by-m numeric array) Array of generated points where
    %            each column represents a point's (x, y) coordinates.
    %            The output size is (2, n*numPoints) where `n` is the 
    %            number of input points.
    %
    % Example:
    %   % Define center points
    %   xy = [0 0; 5 5];
    %   % Define radius and number of points
    %   r = 3;
    %   numPoints = 10;
    %   % Get points within the radius
    %   points = getPointsWithinRadius(xy, r, numPoints);
    %
    % Notes:
    %   - The function removes duplicate points from the output (optional).
    %   - If the input `xy` does not have two columns, an error is thrown.
    %
    % See also: linspace, unique, cos, sin
    
    % Validate input dimensions
    if size(xy, 2) ~= 2
        error('xy coordinates array shape should be n-by-2.');
    end
    
    % Initialize an empty array to hold the points
    points = [];
    
    % Loop through each point defined by (x_center, y_center)
    for i = 1:size(xy, 1)
        % Current point coordinates
        x_center = xy(i, 1);
        y_center = xy(i, 2);
        
        % Generate points around the current center point
        theta = linspace(0, 2 * pi, numPoints); % Generate numPoints around the circle
        x_circle = x_center + r * cos(theta);
        y_circle = y_center + r * sin(theta);
        
        % Combine x and y coordinates into a 2-by-m array
        new_points = [x_circle; y_circle];
        
        % Append to the results
        points = [points, new_points]; %#ok<AGROW> 
    end
    
    % Remove duplicate points (optional)
    points = unique(points', 'rows')'; % Transpose for unique, then back
end

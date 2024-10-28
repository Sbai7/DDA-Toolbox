function [hfig, hbar] = plotInjectionCatchmentPartitions(well_pair, well_pair_volume, colors, inj_wells, pump_wells)
% Inputs:
% well_pair - Nx2 cell array of well pairs (e.g., {'I1', 'P3'; 'I1', 'P4'; ...})
% well_pair_volume - Nx1 array of volumes corresponding to each well pair
% colors - Nx3 matrix of RGB values for each well pair
% inj_wells - Cell array of injection well names (e.g., {'I1', 'I2', ...})
% pump_wells - Cell array of pumping well names (e.g., {'P1', 'P2', ...})

% Get the number of injection wells
N_inj = numel(inj_wells);
    
% Initialize data structures
total_volumes = zeros(N_inj, 1);   % Total volume for each injection well
contribution_matrix = zeros(N_inj, numel(pump_wells));  % Contributions from each pumping well
color_matrix = zeros(N_inj, numel(pump_wells));

% Calculate total volumes and contributions
for i = 1:size(well_pair, 1)
   inj_index = find(strcmp(inj_wells, well_pair{i, 1}));
   pump_index = find(strcmp(pump_wells, well_pair{i, 2}));

   if ~isempty(inj_index) && ~isempty(pump_index)
      total_volumes(inj_index) = total_volumes(inj_index) + well_pair_volume(i);
      contribution_matrix(inj_index, pump_index) = contribution_matrix(inj_index, pump_index) + well_pair_volume(i);
      color_matrix(inj_index, pump_index) = i;
   end
end

% Pack the color matrix 
color_matrix = color_matrix(:);

M = length(color_matrix);
newColors = ones(M, 3); % Pre-allocate M-by-3 with white color
    
% Assign colors from the colors array for non-zero indices
for i = 1:M
   if color_matrix(i) ~= 0
      newColors(i, :) = colors(color_matrix(i), :);  % Get RGB color from colors
   end
end

% Create the horizontal stacked bar plot
hfig = figure;
hbar = barh(1:1:N_inj, contribution_matrix, 'stacked', 'FaceColor', 'flat');
hold on;

% Assign custom colors to each segment of the bars
% We have 3 segments, and each segment has 2 categories (1980 and 1990)
for k = 1:length(hbar)
   % Each bar will take its color from the respective segments
   % Since h(k) corresponds to the kth segment across all categories, we
   % assign colors for each segment for all bars
   hbar(k).CData = newColors((k-1)*N_inj + (1:N_inj), :);
end

% Add the values inside each bar segment
for i = 1:N_inj  % Loop through each bar (category)
   cumulativeWidth = 0;  % Track cumulative width for each segment
    
   for j = 1:size(contribution_matrix, 1)  % Loop through each segment in the bar
      % Get the width of the current segment (absolute value for correct positioning)
      segmentWidth = contribution_matrix(i,j);
      value = 100 * contribution_matrix(i,j) / sum(contribution_matrix(i,:));
        
      % Center the text in the segment (x-coordinate)
      textPosX = cumulativeWidth + segmentWidth / 2;
        
        % Y-coordinate is the category value
        textPosY = i;
        
        % Display the text (use num2str to convert value to string)
        if segmentWidth ~= 0
            text(textPosX, textPosY, [num2str(value, 2) '%'], ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'Color', 'w');
        end
        
        % Update cumulative width for the next segment
        cumulativeWidth = cumulativeWidth + segmentWidth;
    end
end


% Customize the plot
set(gca, 'YTick', 1:N_inj, 'YTickLabel', inj_wells);
grid on;
box on;
hold off;

end

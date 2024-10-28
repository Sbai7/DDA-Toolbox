function [sorted_rt, cumulative_pore_volume, cumulative_fractional_q] = ...
    computeCumulativePhiQ(pore_volume, residence_time)
%computeCumulativePhiQ: Computes cumulative pore volumes and fractional
%flow rates along increasing values of a residence time distribution. 
%
% SYNOPSIS:
%   [sorted_rt, cumulative_pore_volume, cumulative_fractional_q] = ...
%          computeCumulativePhiQ(pore_volume, residence_time)
%
% DESCRIPTION:
%
% PARAMETERS:
%   pore_volume    - pore volume, one value per cell
%   residence_time - residence time array, one value per cell
%
% RETURNS:
%   cumulative_pore_volume - cumulative pore volumes for increasing
%   residence time values
%
%   cumulative_fractional_q - cumulative flow rate for increasing residence
%   time values.
%

assert(size(residence_time,2) == 1, 'Residence time array must have one column.')

% Sort cells based on residence time.
[sorted_rt,order] = sort(residence_time);  

% Sort pore volumes in same order of residence time values.
sorted_pore_volume = pore_volume(order);

% Compute cumulative sum of sorted pore volumes.
cumulative_pore_volume = cumsum(sorted_pore_volume);

% Normalize the cumulative pore volumes array.
total_volume = full(cumulative_pore_volume(end));
cumulative_pore_volume = cumulative_pore_volume/total_volume;

% Compute cumulative sum of fractional flow rates.
sorted_fractional_q = sorted_pore_volume./sorted_rt;
cumulative_fractional_q = cumsum(sorted_fractional_q);
total_q = full(cumulative_fractional_q(end));

% Normalize the fractional flow rates array.
cumulative_fractional_q = cumulative_fractional_q/total_q; 

end
function [r,SF3,S3L_radial,S3T_radial, ...
    radial_pair_count,radial_direction_count] = ...
    calc_radial_nonperiodic( ...
        S3L_map,S3T_map,pair_count, ...
        dx_mean,dy_mean,r_edges)

% CALC_RADIAL_NONPERIODIC
% Radially averages nonperiodic third-order structure functions
% using prescribed separation-bin boundaries.
%
% S3L_map and S3T_map are assumed to be normalized by pair_count
% at each displacement vector before being passed to this function.

[ny_lag,nx_lag] = size(S3L_map);

if mod(nx_lag,2) ~= 1 || mod(ny_lag,2) ~= 1
    error('The displacement maps must have odd dimensions.');
end

% Physical displacement coordinates
lag_x = (-floor(nx_lag/2):floor(nx_lag/2)) .* dx_mean;
lag_y = (-floor(ny_lag/2):floor(ny_lag/2)) .* dy_mean;

[X,Y] = meshgrid(lag_x,lag_y);
R = sqrt(X.^2 + Y.^2);

n_bins = numel(r_edges)-1;

% Representative radius of each logarithmic bin
r = sqrt(r_edges(1:end-1) .* r_edges(2:end));

S3L_radial = nan(1,n_bins);
S3T_radial = nan(1,n_bins);
SF3 = nan(1,n_bins);

radial_pair_count = zeros(1,n_bins);
radial_direction_count = zeros(1,n_bins);

for i = 1:n_bins

    if i < n_bins
        radial_mask = ...
            R >= r_edges(i) & R < r_edges(i+1);
    else
        radial_mask = ...
            R >= r_edges(i) & R <= r_edges(i+1);
    end

    valid = radial_mask & ...
            isfinite(S3L_map) & ...
            isfinite(S3T_map) & ...
            isfinite(pair_count) & ...
            pair_count > 0;

    radial_direction_count(i) = nnz(valid);

    if radial_direction_count(i) == 0
        continue
    end

    weights = double(pair_count(valid));
    total_pairs = sum(weights);

    radial_pair_count(i) = total_pairs;

    if total_pairs <= 0
        continue
    end

    % Pair-count-weighted radial averages
    S3L_radial(i) = ...
        sum(double(S3L_map(valid)) .* weights) ./ total_pairs;

    S3T_radial(i) = ...
        sum(double(S3T_map(valid)) .* weights) ./ total_pairs;

    SF3(i) = S3L_radial(i) + S3T_radial(i);

end

end
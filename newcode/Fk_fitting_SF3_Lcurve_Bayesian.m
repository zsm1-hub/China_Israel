function [SpecFlux, Vt, ebs, kf, lf, lambda_opt] = Fk_fitting_SF3_Lcurve_Bayesian(SF3, ...
    dist_axis, mindist, maxdist, kftype, inv_style, lambda_vec, kf1,plt)
% Fk_fitting_SF3_Lcurve - 
%
% input:
%   SF3:  (r × time)
%   dist_axis:  (1 × r)
%   mindist:  choose smallest range (m)
%   maxdist:  (m)
%   kftype:  ('log' or 'linear')
%   inv_style:  ('NNLS', 'LS', 'RLS')
%   lambda_vec: Regularized vec (for RLS)
%   kf1: didn't recommend (if you don't use can make kf1=1)
%
% output:
%   SpecFlux: fk (m²/s³)
%   Vt: fitting sf3
%   ebs: energy injection rate (m³/s³)
%   kf:  (m⁻¹)
%   lf:  (m)
%   lambda_opt: 

% 
% nsamps = size(SF3, 2);
r = dist_axis;
Nr = length(r);

% 
ns = find(dist_axis >= mindist, 1);
ne = find(dist_axis <= maxdist, 1, 'last');

if isempty(ns) || isempty(ne)
    error('No valid points found in the specified distance range.');
end

R = r(ns:ne);
NR = length(R);
lf = R;

% 
if strcmp(kftype, 'log')
    kf = logspace(log10(1/max(R)), log10(1/min(R)), length(R)-1).*2.*pi;
else
    kf = kf1;
end

dk = diff(kf);
kf = 0.5*(kf(1:end-1) + kf(2:end));
kf = fliplr(kf);
dk = fliplr(dk);

Nk = length(kf);

% 
ebs = zeros(Nk+1, 1);
S = zeros(Nr, 1);
Vt = zeros(NR, 1);
SpecFlux = zeros(Nk, 1);




% 
S(:, 1) = SF3(:, 1);
V = S(ns:ne, 1)';

% make matrix Ax=b
A = zeros(NR, Nk+1);
for j = 1:Nk
    A(:, j) = -4/kf(j) * besselj(1, kf(j)*R)' * dk(j);
end
A(:, end) = 2*R';

% 
norm_flag = 1;
if norm_flag == 1
    for j = 1:NR
        A(j, :) = A(j, :) ./ abs(R(j));
    end
    V = V ./ abs(R);
end

% 
if strcmp(inv_style, 'NNLS')
    ebs(:, 1) = lsqnonneg(A, V');
elseif strcmp(inv_style, 'LS')
    ebs(:, 1) = A \ V';
elseif strcmp(inv_style, 'RLS')
    % L curve

    lambda_opt = lambda_vec(1);
end
% 
n_cols = size(A, 2);
A_aug = [A; sqrt(lambda_opt) * eye(n_cols)];
V_aug = [V'; zeros(n_cols, 1)];
ebs(:, 1) = A_aug \ V_aug;

% 
Vt(:, 1) = 2 * ebs(end, 1) * R;
for j = 1:Nk
    Vt(:, 1) = Vt(:, 1) - 4 * ebs(j, 1) ./ kf(j) .* besselj(1, kf(j)*R)' * dk(j);
end

% 
SpecFlux(1, 1) = -ebs(end, 1);
for j = 2:Nk
    SpecFlux(j, 1) = SpecFlux(j-1, 1) + ebs(end-j+1, 1) * dk(end-j+2);
end


% 
SpecFlux = flipud(SpecFlux);
return
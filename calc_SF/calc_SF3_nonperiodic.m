function [S3L, S3T, pair_count, X, Y] = calc_SF3_nonperiodic(u, v)
%CALC_SF3_NONPERIODIC Third-order structure functions without periodic wrap.
%
%   [S3L,S3T,PAIR_COUNT,X,Y] = CALC_SF3_NONPERIODIC(U,V)
%
% U and V are collocated two-dimensional velocity fields. The calculation
% uses zero-padded FFT correlations, retains all one-point terms in the
% velocity-increment expansions, and normalizes every displacement by its
% actual number of valid point pairs. NaNs in U or V are treated as invalid
% grid points. No periodicity or statistical homogeneity is invoked in the
% calculation of the direction-dependent moments.
%
% Correlation convention:
%   C_ab(r) = < a(x+r) b(x) >_valid pairs.
%
% S3L is <delta u_L^3>; S3T is <delta u_L delta u_T^2>.

    if ~isequal(size(u), size(v))
        error('u and v must have the same size.');
    end
    if ~ismatrix(u) || isempty(u)
        error('u and v must be nonempty two-dimensional arrays.');
    end

    u = double(u);
    v = double(v);
    valid = isfinite(u) & isfinite(v);

    u0 = zeros(size(u));
    v0 = zeros(size(v));
    u0(valid) = u(valid);
    v0(valid) = v(valid);
    one = double(valid);

    [ny, nx] = size(u0);
    py = 2*ny - 1;
    px = 2*nx - 1;

    % Fourier transforms needed by the complete cubic expansions.
    F1   = fft2(one,       py, px);
    Fu   = fft2(u0,        py, px);
    Fv   = fft2(v0,        py, px);
    Fu2  = fft2(u0.^2,     py, px);
    Fv2  = fft2(v0.^2,     py, px);
    Fuv  = fft2(u0.*v0,    py, px);
    Fu3  = fft2(u0.^3,     py, px);
    Fv3  = fft2(v0.^3,     py, px);
    Fvu2 = fft2(v0.*u0.^2, py, px);
    Fuv2 = fft2(u0.*v0.^2, py, px);

    pair_count = real(fftshift(ifft2(F1 .* conj(F1))));
    % Remove roundoff around integer pair counts.
    pair_count = round(max(pair_count, 0));

    % C(Fa,Fb) represents <a(x+r)b(x)> over valid pairs.
    C = @(Fa,Fb) local_corr(Fa, Fb, pair_count);

    % Complete expansions: no cancellation based on homogeneity is used.
    Suuu = C(Fu3,F1) - 3*C(Fu2,Fu) + 3*C(Fu,Fu2) - C(F1,Fu3);
    Svvv = C(Fv3,F1) - 3*C(Fv2,Fv) + 3*C(Fv,Fv2) - C(F1,Fv3);

    % Svuu = <delta v (delta u)^2>.
    Svuu = C(Fvu2,F1) - 2*C(Fuv,Fu) + C(Fv,Fu2) ...
          - C(Fu2,Fv) + 2*C(Fu,Fuv) - C(F1,Fvu2);

    % Suvv = <delta u (delta v)^2>.
    Suvv = C(Fuv2,F1) - 2*C(Fuv,Fv) + C(Fu,Fv2) ...
          - C(Fv2,Fu) + 2*C(Fv,Fuv) - C(F1,Fuv2);

    [X, Y] = meshgrid(-(nx-1):(nx-1), -(ny-1):(ny-1));
    R = hypot(X, Y);
    costheta = zeros(size(R));
    sintheta = zeros(size(R));
    nonzero = R > 0;
    costheta(nonzero) = X(nonzero)./R(nonzero);
    sintheta(nonzero) = Y(nonzero)./R(nonzero);

    c = costheta;
    s = sintheta;
    S3L = c.^3.*Suuu + 3*c.^2.*s.*Svuu ...
        + 3*c.*s.^2.*Suvv + s.^3.*Svvv;

    S3T = s.^2.*c.*Suuu + (s.^3 - 2*s.*c.^2).*Svuu ...
        + (c.^3 - 2*c.*s.^2).*Suvv + s.*c.^2.*Svvv;

    % The direction of the zero-separation vector is undefined.
    S3L(~nonzero | pair_count < 1) = NaN;
    S3T(~nonzero | pair_count < 1) = NaN;
end


function out = local_corr(Fa, Fb, pair_count)
% Linear (noncircular) cross-correlation normalized by valid-pair count.
    numerator = real(fftshift(ifft2(Fa .* conj(Fb))));
    out = NaN(size(numerator));
    good = pair_count >= 1;
    out(good) = numerator(good)./pair_count(good);
end

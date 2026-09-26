function r = eq_fourier_eval(theta, B, rmin, dth_half)
%EQ_FOURIER_EVAL  The paper's Fourier-series radius, eq. (40).
%
%   r = EQ_FOURIER_EVAL(theta, B, rmin, dth_half)
%
%       r(theta) = r_min / ( 1 + sum_i B(i+1) * cos(i*pi*theta/dth_half) )
%
%   theta    : polar angle measured from periapsis (any size array)
%   B        : [B0, B1, ...]  (B0 is the constant term)
%   dth_half : angle from periapsis to apoapsis (half of dtheta)
%
%   Order 1:  B = [B0 B1]      (paper eq. 33/36)
%   Order 2:  B = [B0 B1 B2]   (paper Table 1)

phi = pi*theta/dth_half;
D   = ones(size(theta));
for i = 0:numel(B) - 1
    D = D + B(i+1)*cos(i*phi);
end
r = rmin ./ D;
end

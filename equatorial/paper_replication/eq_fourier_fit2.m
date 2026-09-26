function [B, J, rapp] = eq_fourier_fit2(theta, r, rmin, rmax, dth_half)
%EQ_FOURIER_FIT2  Order-2 Fourier fit of r(theta), as in Wang et al. (2014), Table 1.
%
%   [B, J, rapp] = EQ_FOURIER_FIT2(theta, r, rmin, rmax, dth_half)
%
%   Model (paper eq. 40, m = 2):
%       r_app = r_min / (1 + B0 + B1 cos(phi) + B2 cos(2 phi)),  phi = pi*theta/dth_half
%
%   Constraints (paper eq. 41): r_app equals r_min at theta = 0 and r_max at
%   theta = dth_half. Two equations in three unknowns, which gives
%       B1 = -q/2,   B0 = q/2 - B2,   with q = r_min/r_max - 1
%   so only B2 is free. It is chosen to minimise the paper's objective (42):
%       J = sum over samples of (r_app - r)^2
%
%   Because only one number is fitted, a few Gauss-Newton steps (Newton's
%   method for least squares) reach the minimum.
%
%   Inputs : theta, r = exact samples over the half period (from integration)
%   Outputs: B = [B0 B1 B2], J = objective value, rapp = fitted radius

q   = rmin/rmax - 1;
B1  = -q/2;
phi = pi*theta(:)/dth_half;
rr  = r(:);

B2 = 0;
for k = 1:50
    D    = 1 + q/2 + B1*cos(phi) + B2*(cos(2*phi) - 1);
    res  = rmin./D - rr;
    dr   = -rmin*(cos(2*phi) - 1)./D.^2;     % d r_app / d B2
    step = -(dr.'*res)/(dr.'*dr);
    B2   = B2 + step;
    if abs(step) < 1e-18
        break
    end
end

B    = [q/2 - B2, B1, B2];
rapp = eq_fourier_eval(theta, B, rmin, dth_half);
J    = sum((rapp(:) - rr).^2);
end

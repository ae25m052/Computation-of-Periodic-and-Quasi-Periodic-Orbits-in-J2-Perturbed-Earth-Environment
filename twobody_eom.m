function dxdt = twobody_eom(~, x, mu, epsilon, J2, Re)
%   dxdt = twobody_eom(t, x, mu, epsilon, J2, Re)
%   Inputs:
%       x        - state vector [r(3); v(3)]  (non-dimensional if mu is)
%       mu       - gravitational parameter (use 1 if non-dimensionalized)
%       epsilon  - continuation parameter, 0 (pure Keplerian) -> 1 (full J2)
%       J2       - Earth's J2 coefficient (dimensionless, ~1.08263e-3)
%       Re       - equatorial radius, SAME UNITS as r in x (non-dim: Re/DU)
%
%   Output:
%       dxdt     - time derivative of state vector (6x1)

    r_vec = x(1:3);
    v_vec = x(4:6);
    r = norm(r_vec);

    % --- Point-mass (Keplerian) acceleration ---
    a_pm = -mu / r^3 * r_vec;

    % --- J2 perturbation acceleration (zero when epsilon = 0) ---
    if epsilon ~= 0
        x_ = r_vec(1); y_ = r_vec(2); z_ = r_vec(3);
        factor = -1.5 * J2 * mu * Re^2 / r^5;
        zr2 = (z_/r)^2;
        a_J2 = factor * [ x_ * (1 - 5*zr2);
                           y_ * (1 - 5*zr2);
                           z_ * (3 - 5*zr2) ];
    else
        a_J2 = zeros(3,1);
    end

    a_total = a_pm + epsilon * a_J2;

    dxdt = [v_vec; a_total];
end

function dXdt = stm_eom(t, X, mu, epsilon, J2, Re)
%STM_EOM  Combined state and state-transition-matrix equations of motion.
%
%   dXdt = stm_eom(t, X, mu, epsilon, J2, Re)
%
%   Integrates the 6 state equations and the 36 variational equations
%   simultaneously, i.e. 42 coupled ODEs:
%
%       xdot   = f(x)                        (6  equations)
%       Phidot = A(t) Phi ,  Phi(0) = I      (36 equations)
%
%   where A(t) = [ 0  I ; G  0 ] and G = da/dr is supplied by
%   gravity_gradient.m, evaluated ALONG the trajectory.
%
%   Input state layout:
%       X(1:6)   - [r; v]
%       X(7:42)  - Phi, stored COLUMN-WISE (MATLAB reshape order)
%
%   WHY INTEGRATE TOGETHER RATHER THAN AFTERWARDS
%   A(t) depends on the position at time t, so the variational equations
%   cannot be integrated without the trajectory. Propagating both in one
%   call guarantees the STM sees exactly the same discrete trajectory the
%   integrator produced, rather than an interpolation of it. Interpolating
%   the state and integrating Phi separately introduces an error that is
%   invisible over one period but corrupts the monodromy matrix and hence
%   the Floquet multipliers in Stage 5.
%
%   Initialise with:
%       X0 = [x0; reshape(eye(6), 36, 1)];

    r_vec = X(1:3);
    v_vec = X(4:6);
    Phi   = reshape(X(7:42), 6, 6);

    % ---- state derivative (identical physics to twobody_eom.m) ----
    r = norm(r_vec);
    a_vec = -mu/r^3 * r_vec;

    if epsilon ~= 0
        x_ = r_vec(1); y_ = r_vec(2); z_ = r_vec(3);
        factor = -1.5 * J2 * mu * Re^2 / r^5;
        zr2 = (z_/r)^2;
        a_J2 = factor * [ x_ * (1 - 5*zr2);
                          y_ * (1 - 5*zr2);
                          z_ * (3 - 5*zr2) ];
        a_vec = a_vec + epsilon * a_J2;
    end

    % ---- variational equations ----
    G = gravity_gradient(r_vec, mu, epsilon, J2, Re);

    A = [ zeros(3), eye(3);
          G,        zeros(3) ];

    Phidot = A * Phi;

    dXdt = [ v_vec;
             a_vec;
             reshape(Phidot, 36, 1) ];
end

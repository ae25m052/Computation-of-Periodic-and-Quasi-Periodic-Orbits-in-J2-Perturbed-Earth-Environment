function [XT, Phi] = prop_stm(X0, tf, epsilon, tol)
%PROP_STM  Propagate state and state transition matrix (the 42 equations).
%
%   [XT, Phi] = prop_stm(X0, tf, epsilon, tol)
%
%   Inputs:
%       X0       initial state [r; v], canonical, 6x1
%       tf       final time, canonical -- pass the KEPLERIAN period
%       epsilon  0 = pure Kepler, 1 = full J2
%       tol      RelTol and AbsTol for ode89
%
%   Outputs:
%       XT       state at tf, 6x1
%       Phi      state transition matrix at tf, 6x6
%
%   Integrates
%       rdot   = v
%       vdot   = -mu*r/|r|^3 + epsilon * (J2 acceleration)
%       Phidot = A(r, epsilon) * Phi ,    Phi(0) = I
%
%   through the project's j2_stm_eom.m, which builds A at every step from
%   the analytically derived j2_A_matrix.m -- the nine hand-derived entries
%   of dG/dr.
%
%   IMPORTANT -- WHAT Phi(T) IS AND IS NOT
%   At epsilon = 1 the orbit is NOT periodic, so Phi(T) is the STM over
%   that time interval and is NOT a monodromy matrix. Its eigenvalues are
%   therefore NOT Floquet multipliers and do not measure orbital stability.
%   A large value simply reflects a defective (Jordan-block) Keplerian STM
%   responding to a perturbation; defective eigenvalues split like sqrt of
%   the perturbation rather than linearly. Floquet multipliers only become
%   meaningful once a genuine periodic orbit has been computed and Phi is
%   taken over that orbit's own period.
%
%   M.Tech. project, Raj Khismatrao (AE25M052).

    c    = j2_constants();
    opts = odeset('RelTol', tol, 'AbsTol', tol);

    Y0 = [X0(:); reshape(eye(6), 36, 1)];

    [~, Y] = ode89(@(t,Y) j2_stm_eom(t, Y, c.mu, epsilon, c.J2, c.Re), ...
                   [0 tf], Y0, opts);

    XT  = Y(end, 1:6).';
    Phi = reshape(Y(end, 7:42).', 6, 6);
end

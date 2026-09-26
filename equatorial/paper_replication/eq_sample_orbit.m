function [t, theta, r, X] = eq_sample_orbit(X0, T, n, mu, J2, Re, opts)
%EQ_SAMPLE_ORBIT  Propagate an orbit and return n evenly spaced samples in time.
%
%   [t, theta, r, X] = EQ_SAMPLE_ORBIT(X0, T, n, mu, J2, Re, opts)
%
%   t     : n x 1 times from 0 to T
%   theta : n x 1 polar angle in the x-y plane, UNWRAPPED (keeps counting
%           past 2*pi instead of jumping back to 0)
%   r     : n x 1 distance from Earth's centre
%   X     : n x 6 states
%
%   Uses j2_stm_eom.m with eps = 1 (the STM part is carried along but not
%   returned).

rhs    = @(tt, Y) j2_stm_eom(tt, Y, mu, 1, J2, Re);
Y0     = [X0(:); reshape(eye(6), 36, 1)];
tspan  = linspace(0, T, n);
[t, Y] = ode89(rhs, tspan, Y0, opts);
X      = Y(:, 1:6);
theta  = unwrap(atan2(X(:,2), X(:,1)));
r      = sqrt(sum(X(:,1:3).^2, 2));
end

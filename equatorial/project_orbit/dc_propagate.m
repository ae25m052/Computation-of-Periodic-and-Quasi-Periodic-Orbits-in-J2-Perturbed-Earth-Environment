function [t, X, PhiEnd] = dc_propagate(X0, tspan, mu, J2, Re, opts)
%DC_PROPAGATE  Integrate state + STM with j2_stm_eom.m (eps = 1).
%
%   [t, X, PhiEnd] = DC_PROPAGATE(X0, tspan, mu, J2, Re, opts)
%
%   tspan  : [t0 tf], or a vector of times at which to return the state
%   t      : output times (column)
%   X      : states, one row per time (N x 6)
%   PhiEnd : the 6x6 STM at the last time

Y0 = [X0(:); reshape(eye(6), 36, 1)];
[t, Y] = ode89(@(tt, Y) j2_stm_eom(tt, Y, mu, 1, J2, Re), tspan, Y0, opts);
X      = Y(:, 1:6);
PhiEnd = reshape(Y(end, 7:42), 6, 6);
end

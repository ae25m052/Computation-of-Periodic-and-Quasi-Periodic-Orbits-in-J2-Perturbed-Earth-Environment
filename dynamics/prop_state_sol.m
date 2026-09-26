function sol = prop_state_sol(X0, tf, epsilon, tol)
%PROP_STATE_SOL  Propagate the 6 state equations, returning a SOLUTION STRUCTURE.
%
%   sol = prop_state_sol(X0, tf, epsilon, tol)
%
%   Inputs:
%       X0       initial state [r; v], canonical, 6x1
%       tf       final time, canonical
%       epsilon  0 = pure Kepler, 1 = full J2
%       tol      RelTol and AbsTol for ode89
%
%   Output:
%       sol      MATLAB ODE solution structure
%
%   WHY A SOLUTION STRUCTURE RATHER THAN SAMPLED ARRAYS
%   Returning sol means the trajectory is integrated ONCE and can then be
%   sampled at any set of times afterwards with deval:
%
%       Xrev = deval(sol, (0:N)*T);     % once per revolution, for elements
%       sol.x, sol.y                    % the integrator's own steps, for drawing
%
%   It also lets ode89 choose its own step sizes, which clusters them near
%   perigee where the motion is fastest. Forcing a uniform time grid would
%   under-sample exactly the part of an e = 0.65 orbit that matters most.
%
%   Uses the project's twobody_eom.m, so the physics here is identical to
%   every earlier stage.
%
%   M.Tech. project, Raj Khismatrao (AE25M052).

    c    = j2_constants();
    opts = odeset('RelTol', tol, 'AbsTol', tol);

    sol = ode89(@(t,x) twobody_eom(t, x, c.mu, epsilon, c.J2, c.Re), ...
                [0 tf], X0(:), opts);
end

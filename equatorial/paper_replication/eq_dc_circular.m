function [r0, T, hist] = eq_dc_circular(h, r0, mu, J2, Re, opts)
%EQ_DC_CIRCULAR  Differential correction for the J2 circular equatorial orbit.
%
%   [r0, T, hist] = EQ_DC_CIRCULAR(h, r0_guess, mu, J2, Re, opts)
%
%   Finds the radius r0 of the circular orbit in the equatorial plane
%   (i = 0) for a FIXED angular momentum h. This is the paper's "critical
%   circular orbit" (Wang et al. 2014, Sec. 4.2). It is a genuinely
%   periodic orbit in ECI: after one period T the state repeats exactly.
%
%   ---- The idea (symmetric single shooting) ---------------------------
%   Start on the +x axis with the velocity pointing straight along +y:
%       X0 = [r0; 0; 0;  0; h/r0; 0]
%   (xdot0 = 0 means we start at an apsis, i.e. rdot = 0.)
%   Integrate until the orbit crosses the -x axis (half a revolution).
%   If the orbit is circular it crosses the -x axis perpendicularly, so
%       g(r0) = xdot at the crossing = 0.
%   A non-circular orbit reaches its opposite apsis slightly AFTER
%   180 deg (the J2 apsidal advance), so it crosses the axis with
%   xdot ~= 0. That is the error we drive to zero.
%
%   ---- The Newton step -------------------------------------------------
%   dX0/dr0 = [1; 0; 0; 0; -h/r0^2; 0]   (h is held fixed)
%   A change in X0 also changes WHEN the crossing happens. Correcting for
%   that (the crossing is where y = 0) gives
%       dg/dX0 = Phi(4,:) - (xddot_f / ydot_f) * Phi(2,:)
%   and then  r0_new = r0 - g / (dg/dX0 * dX0/dr0).
%
%   Inputs
%     h        : angular momentum (canonical), e.g. 1.1205
%     r0       : starting guess for the radius (canonical), e.g. h^2
%     mu, J2, Re : constants (mu = 1, Re = 1 in canonical units)
%     opts     : odeset options (tolerances). Events are added here.
%
%   Outputs
%     r0   : corrected radius
%     T    : full period = 2 x (time to the half-revolution crossing)
%     hist : one row per iteration: [iteration, r0, g, dg/dr0]
%
%   Calls j2_stm_eom.m (the validated 42-equation right-hand side) with
%   eps = 1, and eq_event_ycross.m.

tol    = 1e-12;       % |xdot| at the crossing, canonical velocity units
it_max = 20;
optsE  = odeset(opts, 'Events', @eq_event_ycross);
hist   = zeros(0, 4);
rhs    = @(t, Y) j2_stm_eom(t, Y, mu, 1, J2, Re);

for k = 1:it_max
    X0   = [r0; 0; 0; 0; h/r0; 0];
    Y0   = [X0; reshape(eye(6), 36, 1)];
    tmax = 4*pi*sqrt(r0^3/mu);                 % two Kepler periods: ample
    [~, ~, te, ye] = ode89(rhs, [0 tmax], Y0, optsE);
    if isempty(te)
        error('eq_dc_circular:noCrossing', 'No x-axis crossing found.');
    end
    tf  = te(end);
    Yf  = ye(end, :).';
    Xf  = Yf(1:6);
    Phi = reshape(Yf(7:42), 6, 6);

    dYf = j2_stm_eom(tf, Yf, mu, 1, J2, Re);   % gives the acceleration
    af  = dYf(4:6);

    g   = Xf(4);                               % xdot at the crossing
    dX0 = [1; 0; 0; 0; -h/r0^2; 0];
    dg  = (Phi(4,:) - af(1)/Xf(5)*Phi(2,:)) * dX0;

    hist(end+1, :) = [k-1, r0, g, dg]; %#ok<AGROW>
    if abs(g) < tol
        break
    end
    r0 = r0 - g/dg;
end
if abs(g) >= tol
    warning('eq_dc_circular:notConverged', ...
        'Not converged after %d iterations, |g| = %.2e', it_max, abs(g));
end
T = 2*tf;
end

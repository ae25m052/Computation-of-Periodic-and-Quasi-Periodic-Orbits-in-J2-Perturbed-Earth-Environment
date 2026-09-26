function out = eq_dc_pseudo(h, e, mu, J2, Re, opts)
%EQ_DC_PSEUDO  Differential correction for a J2 "pseudo-elliptic" equatorial orbit.
%
%   out = EQ_DC_PSEUDO(h, e, mu, J2, Re, opts)
%
%   For a FIXED angular momentum h and a target J2-eccentricity e, finds the
%   periapsis radius r_min so that the orbit's apoapsis is exactly
%       r_max = r_min * (1 + e)/(1 - e).
%   This is the paper's definition e_J2 = (r_max - r_min)/(r_max + r_min),
%   eq. (35) of Wang et al. (2014).
%
%   ---- What kind of orbit this is --------------------------------------
%   In ECI this orbit is NOT periodic. Its shape repeats every radial period
%   T_r, but by then the ellipse has turned a little (the apsidal advance).
%   So the exact statement is
%       X(T_r) = R_z(dtheta) * X(0)
%   where R_z is a rotation about the z-axis by the angle dtheta covered in
%   one radial period (slightly more than 2*pi). Periodic in a frame that
%   turns with the ellipse; quasi-periodic (two frequencies) in ECI.
%
%   ---- The shooting problem --------------------------------------------
%   Start at periapsis on the +x axis:  X0 = [r_min; 0; 0; 0; h/r_min; 0].
%   Integrate to apoapsis (eq_event_apoapsis.m: r.v goes from + to -).
%   Error:  g(r_min) = |r_f| - r_min*(1+e)/(1-e).
%   Newton: dg/dr_min = (r_f/|r_f|) . (Phi(1:3,:) * dX0/dr_min) - (1+e)/(1-e)
%   The usual "the stopping time moves too" term is multiplied by rdot at
%   the stop, and rdot = 0 at apoapsis, so it drops out.
%
%   Outputs (struct)
%     rmin, rmax : periapsis and apoapsis radii
%     Tr         : radial period (periapsis to periapsis)
%     dtheta     : polar angle covered in one radial period
%     E          : energy (the paper's epsilon)
%     X0         : corrected initial state
%     hist       : [iteration, r_min, g] per iteration
%
%   Calls j2_stm_eom.m (eps = 1) and eq_event_apoapsis.m.

tol    = 1e-12;
it_max = 30;
ratio  = (1 + e)/(1 - e);

% Starting guess. rc is the J2 circular radius, eq. (21) of the paper;
% rc/(1+e) is the Kepler relation r_min = p/(1+e) with p replaced by rc.
% It is only a guess: the correction below does the real work.
rc   = (h^2 + sqrt(h^4 - 6*J2*mu^2*Re^2))/(2*mu);
rmin = rc/(1 + e);

optsE = odeset(opts, 'Events', @eq_event_apoapsis);
rhs   = @(t, Y) j2_stm_eom(t, Y, mu, 1, J2, Re);
hist  = zeros(0, 3);

for k = 1:it_max
    X0   = [rmin; 0; 0; 0; h/rmin; 0];
    Y0   = [X0; reshape(eye(6), 36, 1)];
    a    = rmin/(1 - e);                          % Kepler-sized estimate
    tmax = 4*pi*sqrt(a^3/mu);                     % two Kepler periods
    [~, ~, te, ye] = ode89(rhs, [0 tmax], Y0, optsE);
    if isempty(te)
        error('eq_dc_pseudo:noApoapsis', 'No apoapsis found.');
    end
    tf  = te(end);
    Yf  = ye(end, :).';
    Xf  = Yf(1:6);
    Phi = reshape(Yf(7:42), 6, 6);

    rf  = norm(Xf(1:3));
    g   = rf - ratio*rmin;
    dX0 = [1; 0; 0; 0; -h/rmin^2; 0];
    dg  = (Xf(1:3).'/rf) * (Phi(1:3,:) * dX0) - ratio;

    hist(end+1, :) = [k-1, rmin, g]; %#ok<AGROW>
    if abs(g) < tol
        break
    end
    rmin = rmin - g/dg;
end
if abs(g) >= tol
    warning('eq_dc_pseudo:notConverged', ...
        'Not converged after %d iterations, |g| = %.2e', it_max, abs(g));
end

th_half = mod(atan2(Xf(2), Xf(1)), 2*pi);     % angle of apoapsis

out.rmin   = rmin;
out.rmax   = rf;
out.Tr     = 2*tf;
out.dtheta = 2*th_half;
out.E      = eq_energy(X0, mu, J2, Re);
out.X0     = X0;
out.hist   = hist;
end

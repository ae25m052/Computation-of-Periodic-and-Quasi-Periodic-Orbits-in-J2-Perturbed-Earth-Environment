function out = dc_correct_apogee(rp, ra_target, vp_guess, mu, J2, Re, opts)
%DC_CORRECT_APOGEE  Differential correction: J2 orbit through a given perigee AND apogee.
%
%   out = DC_CORRECT_APOGEE(rp, ra_target, vp_guess, mu, J2, Re, opts)
%
%   Keeps the perigee radius rp fixed and corrects only the perigee speed
%   vp, so that under the full J2 dynamics the orbit reaches exactly
%   ra_target at apoapsis.
%
%   ---- Why a correction is needed --------------------------------------
%   The Kepler speed sqrt(mu(1+e)/rp) gives apogee ra in the two-body
%   problem. J2 adds extra pull in the equatorial plane, so with that same
%   speed the satellite falls short of ra. We need a slightly larger vp.
%
%   ---- The shooting problem (all canonical) ----------------------------
%   Start at perigee:           X0 = [rp; 0; 0; 0; vp; 0]
%   Integrate state + STM to apoapsis (dc_event_apoapsis.m)
%   Miss:                       g(vp) = |r_f| - ra_target
%   Sensitivity from the STM:   dg/dvp = (r_f/|r_f|) . Phi(1:3, 5)
%      (column 5 of Phi = effect of a change in the initial ydot = vp.
%       The stopping time also moves, but that term is multiplied by rdot
%       at the stop, which is zero at apoapsis, so it drops out.)
%   Newton update:              vp_new = vp - g / (dg/dvp)
%
%   Output struct
%     vp, X0      corrected perigee speed and state
%     Tr          radial period: perigee to perigee = 2 x (perigee to apogee)
%     dtheta      polar angle covered in one radial period (> 2*pi)
%     delta       apsidal advance per radial period = dtheta - 2*pi
%     ra          apoapsis radius reached
%     hist        one row per iteration: [iteration, vp, g (miss), r_apoapsis]
%
%   Calls j2_stm_eom.m (eps = 1) and dc_event_apoapsis.m.

tol    = 1e-11;          % DU  (1e-11 DU = 0.06 mm)
it_max = 20;
optsE  = odeset(opts, 'Events', @dc_event_apoapsis);
rhs    = @(t, Y) j2_stm_eom(t, Y, mu, 1, J2, Re);
a_est  = (rp + ra_target)/2;
tmax   = 3*pi*sqrt(a_est^3/mu);          % 1.5 Kepler periods: apoapsis is at ~0.5

vp   = vp_guess;
hist = zeros(0, 4);
for k = 1:it_max
    X0 = [rp; 0; 0; 0; vp; 0];
    [~, ~, te, ye] = ode89(rhs, [0 tmax], [X0; reshape(eye(6), 36, 1)], optsE);
    if isempty(te)
        error('dc_correct_apogee:noApoapsis', 'No apoapsis found.');
    end
    Yf  = ye(end, :).';
    Xf  = Yf(1:6);
    Phi = reshape(Yf(7:42), 6, 6);
    rf  = norm(Xf(1:3));

    g  = rf - ra_target;
    dg = (Xf(1:3).'/rf) * Phi(1:3, 5);

    hist(end+1, :) = [k-1, vp, g, rf]; %#ok<AGROW>
    if abs(g) < tol
        break
    end
    vp = vp - g/dg;
end
if abs(g) >= tol
    warning('dc_correct_apogee:notConverged', 'Not converged, |g| = %.2e', abs(g));
end

out.vp     = vp;
out.X0     = X0;
out.Tr     = 2*te(end);
out.dtheta = 2*mod(atan2(Xf(2), Xf(1)), 2*pi);
out.delta  = out.dtheta - 2*pi;
out.ra     = rf;
out.hist   = hist;
end

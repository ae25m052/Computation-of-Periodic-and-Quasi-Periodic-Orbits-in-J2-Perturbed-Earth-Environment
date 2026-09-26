function S = brouwer_rates(s, inc_deg)
%BROUWER_RATES  First-order secular rates from Brouwer's theory. CITED, not derived.
%
%   S = brouwer_rates(s, inc_deg)     s from seed_orbit_from_perigee
%
%   These two expressions are standard first-order artificial-satellite
%   theory, quoted from Brouwer (1959). Nothing in this project derives
%   them; they are here purely as an INDEPENDENT YARDSTICK against which
%   the measured element drift can be checked.
%
%       dRAAN/dt = -(3/2) n J2 (Re/p)^2 cos(i)
%       dargp/dt =  (3/4) n J2 (Re/p)^2 (5 cos^2(i) - 1)
%
%   with n the mean motion and p = a(1-e^2) the semi-latus rectum.
%
%   Output struct:
%       dRAAN, dargp                  rad per canonical time unit
%       dRAAN_deg_per_day             deg/day
%       dargp_deg_per_day             deg/day
%
%   HOW TO READ THE COMPARISON
%   Agreement to a fraction of a degree over 400 revolutions confirms the
%   propagation. The residual that remains is the SECOND-ORDER J2^2 effect,
%   which first-order theory omits by construction -- so a small
%   disagreement is expected physics, not an error.
%
%   That residual is exactly why the critical inclination is not a clean
%   zero: at i = 63.434949 deg the first-order dargp vanishes identically,
%   yet the measured drift over 400 revolutions is +0.023 deg, and the true
%   frozen inclination for this orbit sits near 63.443 deg. The same point
%   appears in the Route 6 analysis: the first-order degeneracy resolves
%   only at second order.
%
%   M.Tech. project, Raj Khismatrao (AE25M052).

    c = j2_constants();

    i = deg2rad(inc_deg);
    n = 2*pi / s.T;                 % mean motion, canonical
    p = s.p_km / c.DU;              % semi-latus rectum, canonical

    S.dRAAN = -1.5  * n * c.J2 * (c.Re/p)^2 * cos(i);
    S.dargp =  0.75 * n * c.J2 * (c.Re/p)^2 * (5*cos(i)^2 - 1);

    per_day = rad2deg(1) * 86400 / c.TU;
    S.dRAAN_deg_per_day = S.dRAAN * per_day;
    S.dargp_deg_per_day = S.dargp * per_day;
end

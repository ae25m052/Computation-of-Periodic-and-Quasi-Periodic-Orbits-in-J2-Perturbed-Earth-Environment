function X0 = seed_initial_state(s, inc_deg, RAAN_deg, argp_deg, nu_deg)
%SEED_INITIAL_STATE  Canonical 6x1 state from the seed orbit and three angles.
%
%   X0 = seed_initial_state(s, 45, 0, 0, 0)
%
%   Inputs:
%       s          struct from seed_orbit_from_perigee
%       inc_deg    inclination              [deg]
%       RAAN_deg   right ascension of node  [deg]
%       argp_deg   argument of perigee      [deg]
%       nu_deg     true anomaly at epoch    [deg]
%
%   Output:
%       X0   [r; v] in canonical units, 6x1
%
%   Calls the project's own coe2rv.m. Because coe2rv works in whatever
%   units it is handed, the canonical semi-major axis and mu = 1 go in
%   directly and the result needs no rescaling afterwards.
%
%   M.Tech. project, Raj Khismatrao (AE25M052).

    c = j2_constants();

    [r0, v0] = coe2rv(s.a, s.ecc, deg2rad(inc_deg), deg2rad(RAAN_deg), ...
                      deg2rad(argp_deg), deg2rad(nu_deg), c.mu);

    X0 = [r0(:); v0(:)];
end

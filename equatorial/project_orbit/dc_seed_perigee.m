function seed = dc_seed_perigee(alt_p_km, e, Re_km, mu)
%DC_SEED_PERIGEE  Kepler starting state at perigee for an equatorial orbit.
%
%   seed = DC_SEED_PERIGEE(alt_p_km, e, Re_km, mu)
%
%   Orbital elements (the project orbit):
%       perigee altitude = alt_p_km,  e,  i = 0,  RAAN = 0,  argp = 0,  nu = 0
%
%   NORMALISED FIRST: the altitude is converted to canonical units
%   (DU = Re) on the first line, and everything after is canonical
%   (mu = 1, Re = 1).
%
%   With i = RAAN = argp = nu = 0 the satellite starts at perigee on the
%   +x axis, moving in +y. The usual elements -> state conversion (coe2rv)
%   then reduces to one line:
%       r0 = [r_p; 0; 0],   v0 = [0; v_p; 0],   v_p = sqrt(mu*(1+e)/r_p)
%   This is the KEPLER perigee speed. Under J2 it does not produce the
%   intended apogee; dc_correct_apogee.m fixes that.
%
%   Output struct (all canonical except where a name ends in _km):
%     rp, ra      perigee and apogee radius of the Kepler ellipse [DU]
%     a, e, T     semi-major axis [DU], eccentricity, Kepler period [TU]
%     vp          Kepler perigee speed [DU/TU]
%     X0          6x1 starting state
%     rp_km, ra_km, a_km

rp = (Re_km + alt_p_km)/Re_km;          % canonical: divide by DU = Re
a  = rp/(1 - e);
ra = a*(1 + e);
vp = sqrt(mu*(1 + e)/rp);

seed.rp = rp;   seed.ra = ra;   seed.a = a;   seed.e = e;
seed.T  = 2*pi*sqrt(a^3/mu);
seed.vp = vp;
seed.X0 = [rp; 0; 0; 0; vp; 0];
seed.rp_km = rp*Re_km;  seed.ra_km = ra*Re_km;  seed.a_km = a*Re_km;
end

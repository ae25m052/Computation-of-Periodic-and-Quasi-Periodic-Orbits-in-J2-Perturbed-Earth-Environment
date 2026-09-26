function s = seed_orbit_from_perigee(hp_km, ecc)
%SEED_ORBIT_FROM_PERIGEE  Build the seed orbit from perigee altitude and e.
%
%   s = seed_orbit_from_perigee(300, 0.65)
%
%   The guide specified the PERIGEE ALTITUDE and the eccentricity, not the
%   semi-major axis. So a is not a free choice here -- it follows:
%
%       r_p = Re + h_p ,        a = r_p / (1 - e)
%
%   For h_p = 300 km and e = 0.65 this gives a = 19080.39 km, apogee
%   altitude 25104.5 km, and a Keplerian period of 437.16 min (7.286 h).
%
%   Inputs:
%       hp_km   perigee ALTITUDE above the surface  [km]
%       ecc     eccentricity                        [-]
%
%   Output struct:
%       hp_km  ha_km  rp_km  ra_km  a_km  p_km  ecc      [km, -]
%       a      semi-major axis, canonical               [DU]
%       T      Keplerian period, canonical              [TU]
%       T_min  Keplerian period                         [min]
%
%   NOTE ON T. At epsilon = 1 the orbit is not periodic, so T is not a
%   period of anything -- it is a fixed reference span over which the two
%   models can be compared. The guide's instruction "keep T as Kepler"
%   means exactly this.
%
%   M.Tech. project, Raj Khismatrao (AE25M052).

    c = j2_constants();

    s.hp_km = hp_km;
    s.ecc   = ecc;

    s.rp_km = c.Re_km + hp_km;
    s.a_km  = s.rp_km / (1 - ecc);
    s.ra_km = s.a_km * (1 + ecc);
    s.ha_km = s.ra_km - c.Re_km;
    s.p_km  = s.a_km * (1 - ecc^2);

    s.a     = s.a_km / c.DU;
    s.T     = 2*pi*sqrt(s.a^3 / c.mu);
    s.T_min = s.T * c.TU / 60;
end

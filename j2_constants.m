function c = j2_constants()
%J2_CONSTANTS  Physical constants, canonical units and the critical inclination.
%
%   c = j2_constants()
%
%   Canonical units keep every number near 1, which suits the integrator far
%   better than kilometres and seconds:
%       1 DU = Earth equatorial radius
%       1 TU chosen so that mu = 1
%
%   Fields:
%       mu_km   Earth GM                       [km^3/s^2]
%       Re_km   Earth equatorial radius        [km]
%       J2      J2 coefficient                 [-]
%       DU      length unit                    [km]
%       TU      time unit                      [s]
%       VU      velocity unit                  [km/s]
%       mu      gravitational parameter, canonical  (= 1)
%       Re      equatorial radius, canonical        (= 1)
%       i_crit_deg    critical inclination, prograde   [deg]
%       i_crit_retro  critical inclination, retrograde [deg]
%
%   The critical inclination is where the first-order secular rate of the
%   argument of perigee vanishes, i.e. where 5*cos^2(i) - 1 = 0. It is the
%   one inclination a plain 15-degree grid would step straight over, so it
%   is defined here once and inserted into the sweep explicitly.
%
%   Constants were previously repeated in four separate scripts. This file
%   is the single place they live.
%
%   M.Tech. project, Raj Khismatrao (AE25M052).
%   Guide: Prof. Joel George Manathara, IIT Madras.

    c.mu_km = 398600.4418;              % km^3/s^2
    c.Re_km = 6378.1363;                % km
    c.J2    = 1.08262668e-3;            % dimensionless

    c.DU    = c.Re_km;                  % km
    c.TU    = sqrt(c.DU^3 / c.mu_km);   % s
    c.VU    = c.DU / c.TU;              % km/s

    c.mu    = 1;                        % canonical, by construction
    c.Re    = 1;                        % canonical

    c.i_crit_deg   = acosd(1/sqrt(5));  % 63.434949 deg
    c.i_crit_retro = 180 - c.i_crit_deg;% 116.565051 deg
end

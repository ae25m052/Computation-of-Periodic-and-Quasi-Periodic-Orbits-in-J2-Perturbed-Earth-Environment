function [r_vec, v_vec] = coe2rv(a, e, i, RAAN, argp, nu, mu)
%COE2RV Convert classical orbital elements to ECI position/velocity.
%
%   [r_vec, v_vec] = coe2rv(a, e, i, RAAN, argp, nu, mu)
%
%   Inputs (angles in radians):
%       a     - semi-major axis           [km]
%       e     - eccentricity              [-]
%       i     - inclination                [rad]
%       RAAN  - right ascension of ascending node [rad]
%       argp  - argument of periapsis      [rad]
%       nu    - true anomaly               [rad]
%       mu    - gravitational parameter    [km^3/s^2]
%
%   Outputs:
%       r_vec - position vector in ECI     [km]      (3x1)
%       v_vec - velocity vector in ECI     [km/s]    (3x1)
%

    p = a * (1 - e^2);              % semi-latus rectum
    r_mag = p / (1 + e*cos(nu));

    % Position and velocity in perifocal (PQW) frame
    r_pqw = r_mag * [cos(nu); sin(nu); 0];
    v_pqw = sqrt(mu/p) * [-sin(nu); (e + cos(nu)); 0];

    % Rotation matrices (active rotations, then combined as R3(-RAAN)*R1(-i)*R3(-argp))
    cO = cos(RAAN); sO = sin(RAAN);
    ci = cos(i);    si = sin(i);
    cw = cos(argp); sw = sin(argp);

    R11 = cO*cw - sO*sw*ci;
    R12 = -cO*sw - sO*cw*ci;
    R13 = sO*si;

    R21 = sO*cw + cO*sw*ci;
    R22 = -sO*sw + cO*cw*ci;
    R23 = -cO*si;

    R31 = sw*si;
    R32 = cw*si;
    R33 = ci;

    R_pqw2eci = [R11 R12 R13;
                 R21 R22 R23;
                 R31 R32 R33];

    r_vec = R_pqw2eci * r_pqw;
    v_vec = R_pqw2eci * v_pqw;
end

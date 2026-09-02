function [a, e, i, RAAN, argp, nu] = rv2coe(r_vec, v_vec, mu)
%   [a, e, i, RAAN, argp, nu] = rv2coe(r_vec, v_vec, mu)
%   Inputs:
%       r_vec, v_vec - position/velocity (3x1), any consistent units
%       mu           - gravitational parameter, matching units
%
%   Outputs (angles in radians, wrapped to [0, 2*pi)):
%       a     - semi-major axis
%       e     - eccentricity
%       i     - inclination
%       RAAN  - right ascension of ascending node
%       argp  - argument of periapsis
%       nu    - true anomaly

    r = norm(r_vec);
    v = norm(v_vec);

    h_vec = cross(r_vec, v_vec);
    h = norm(h_vec);

    n_vec = cross([0;0;1], h_vec);   
    n = norm(n_vec);

    e_vec = cross(v_vec, h_vec)/mu - r_vec/r;
    e = norm(e_vec);

    energy = v^2/2 - mu/r;
    a = -mu/(2*energy);

    i = acos(h_vec(3)/h);

    RAAN = atan2(n_vec(2), n_vec(1));
    if RAAN < 0, RAAN = RAAN + 2*pi; end

    argp = acos(dot(n_vec, e_vec)/(n*e));
    if e_vec(3) < 0, argp = 2*pi - argp; end

    nu = acos(dot(e_vec, r_vec)/(e*r));
    if dot(r_vec, v_vec) < 0, nu = 2*pi - nu; end
end

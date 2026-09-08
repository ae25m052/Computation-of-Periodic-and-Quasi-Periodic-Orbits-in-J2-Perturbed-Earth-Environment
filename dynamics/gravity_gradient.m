function G = gravity_gradient(r_vec, mu, epsilon, J2, Re)
%GRAVITY_GRADIENT  3x3 gradient of the gravitational acceleration, da/dr.
%   G = gravity_gradient(r_vec, mu, epsilon, J2, Re)
%       d/dt [r; v] = [ v ; a(r) ]
%       A(t) = [  0    I  ]
%              [  G    0  ]      with  G = da/dr  (this function).
%
%   Inputs (same convention as twobody_eom.m):
%       r_vec   - position vector (3x1), non-dimensional
%       mu      - gravitational parameter (1 in canonical units)
%       epsilon - continuation parameter, 0 (Keplerian) -> 1 (full J2)
%       J2      - J2 coefficient
%       Re      - equatorial radius in the same units as r_vec
%
%   Output:
%       G       - 3x3 matrix da/dr

    x = r_vec(1);  y = r_vec(2);  z = r_vec(3);
    r  = norm(r_vec);

    % ---- point-mass (Keplerian) tidal tensor ----
    G = -mu/r^3 * (eye(3) - 3*(r_vec*r_vec.')/r^2);

    % ---- J2 contribution ----
    if epsilon ~= 0
        k  = -1.5 * J2 * mu * Re^2;
        r5 = r^5;  r7 = r^7;  r9 = r^9;

        GJ = zeros(3,3);

        % row 1:  a_x = k ( x/r^5 - 5 x z^2 / r^7 )
        GJ(1,1) = k*( 1/r5 - 5*x*x/r7 - 5*z*z/r7 + 35*x*x*z*z/r9 );
        GJ(1,2) = k*( -5*x*y/r7 + 35*x*y*z*z/r9 );
        GJ(1,3) = k*( -15*x*z/r7 + 35*x*z*z*z/r9 );

        % row 2:  a_y = k ( y/r^5 - 5 y z^2 / r^7 )   (x <-> y symmetry)
        GJ(2,1) = GJ(1,2);
        GJ(2,2) = k*( 1/r5 - 5*y*y/r7 - 5*z*z/r7 + 35*y*y*z*z/r9 );
        GJ(2,3) = k*( -15*y*z/r7 + 35*y*z*z*z/r9 );

        % row 3:  a_z = k ( 3z/r^5 - 5 z^3 / r^7 )
        GJ(3,1) = GJ(1,3);
        GJ(3,2) = GJ(2,3);
        GJ(3,3) = k*( 3/r5 - 30*z*z/r7 + 35*z^4/r9 );

        G = G + epsilon * GJ;
    end
end

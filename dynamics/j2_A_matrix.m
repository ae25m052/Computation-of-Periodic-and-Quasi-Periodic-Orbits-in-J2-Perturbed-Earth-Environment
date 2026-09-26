function [A, G] = j2_A_matrix(r_vec, mu, eps, J2, Re)
%J2_A_MATRIX  Jacobian A = df/dX for the two-body + J2 problem.
%   [A, G] = J2_A_MATRIX(r_vec, mu, eps, J2, Re)
%   Inputs
%     r_vec : 3x1 position vector [x; y; z]
%     mu    : gravitational parameter GM (1 in canonical units)
%     eps   : continuation parameter. eps = 0 gives pure Kepler,
%             eps = 1 gives the full J2 model. Everything in between is
%             the partially-switched-on problem used in Stage 4.
%     J2    : J2 coefficient (1.08262668e-3 for Earth)
%     Re    : equatorial radius (1 in canonical units)
%   Outputs
%     A : 6x6 Jacobian, A = [0 I ; G 0]
%     G : 3x3 gravity-gradient block, G(i,j) = d(f_i)/d(x_j)
%   The equations of motion this differentiates, with k = (3/2)*eps*J2*mu*Re^2:
%     f_x = -mu*x/r^3 - k*( x/r^5 - 5*x*z^2/r^7 )
%     f_y = -mu*y/r^3 - k*( y/r^5 - 5*y*z^2/r^7 )
%     f_z = -mu*z/r^3 - k*( 3*z/r^5 - 5*z^3/r^7 )

x = r_vec(1);
y = r_vec(2);
z = r_vec(3);

r  = sqrt(x*x + y*y + z*z);
r3 = r^3;
r5 = r^5;
r7 = r^7;
r9 = r^9;

k = 1.5 * eps * J2 * mu * Re^2;     

% ---- the nine entries -------------------------------------------------
G = zeros(3,3);

% dfx/dx
G(1,1) = 3*mu*x*x/r5 - mu/r3 ...
    + k*( -1/r5 + 5*x*x/r7 + 5*z*z/r7 - 35*x*x*z*z/r9 );

% dfx/dy
G(1,2) = 3*mu*x*y/r5 ...
    + k*(  5*x*y/r7 - 35*x*y*z*z/r9 );

% dfx/dz
G(1,3) = 3*mu*x*z/r5 ...
    + k*( 15*x*z/r7 - 35*x*z*z*z/r9 );

% dfy/dx  (equals dfx/dy)
G(2,1) = G(1,2);

% dfy/dy
G(2,2) = 3*mu*y*y/r5 - mu/r3 ...
    + k*( -1/r5 + 5*y*y/r7 + 5*z*z/r7 - 35*y*y*z*z/r9 );

% dfy/dz
G(2,3) = 3*mu*y*z/r5 ...
    + k*( 15*y*z/r7 - 35*y*z*z*z/r9 );

% dfz/dx  
G(3,1) = G(1,3);

% dfz/dy  
G(3,2) = G(2,3);

% dfz/dz
G(3,3) = 3*mu*z*z/r5 - mu/r3 ...
    + k*( -3/r5 + 30*z*z/r7 - 35*z*z*z*z/r9 );

% ---- assemble the 6x6 -------------------------------------------------
A = [zeros(3,3), eye(3);
    G,          zeros(3,3)];

end

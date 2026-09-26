function dY = j2_stm_eom(~, Y, mu, eps, J2, Re)
%J2_STM_EOM  Right-hand side for state + STM, two-body plus J2.
%   dY = J2_STM_EOM(t, Y, mu, eps, J2, Re)
%   Y is 42x1:
%       Y(1:3)   position
%       Y(4:6)   velocity
%       Y(7:42)  the 6x6 STM Phi, unrolled column by column
%   Equations integrated:
%       rdot   = v
%       vdot   = -mu*r/|r|^3 + eps * (J2 acceleration)
%       Phidot = A(r, eps) * Phi,      Phi(0) = I

r_vec = Y(1:3);
v_vec = Y(4:6);
Phi   = reshape(Y(7:42), 6, 6);

x = r_vec(1);
y = r_vec(2);
z = r_vec(3);
r = sqrt(x*x + y*y + z*z);

k = 1.5 * eps * J2 * mu * Re^2;

% ---- state equations ---------------------------------------------------
a_vec = [ -mu*x/r^3 - k*( x/r^5 - 5*x*z^2/r^7 );
    -mu*y/r^3 - k*( y/r^5 - 5*y*z^2/r^7 );
    -mu*z/r^3 - k*( 3*z/r^5 - 5*z^3/r^7 ) ];

% ---- variational equations --------------------------------------------
A    = j2_A_matrix(r_vec, mu, eps, J2, Re);
dPhi = A * Phi;

dY = [v_vec;
    a_vec;
    dPhi(:)];

end

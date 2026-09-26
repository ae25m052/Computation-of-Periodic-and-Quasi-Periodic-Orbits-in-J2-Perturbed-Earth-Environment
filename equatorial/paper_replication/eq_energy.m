function E = eq_energy(X, mu, J2, Re)
%EQ_ENERGY  Specific orbital energy (per unit mass) with the J2 term.
%
%   E = EQ_ENERGY(X, mu, J2, Re)
%
%   X  : 6x1 state [x; y; z; xdot; ydot; zdot]  (canonical units)
%   E  : kinetic + potential energy
%
%   Potential (the same one j2_stm_eom.m differentiates):
%       V = -mu/r + mu*J2*Re^2*(3*z^2/r^2 - 1)/(2*r^3)
%   In the equatorial plane (z = 0) this becomes
%       V = -mu/r - J2*mu*Re^2/(2*r^3)
%   which is equation (2) of Wang et al. (2014).
%
%   In canonical units (mu = 1, Re = 1) E is the paper's epsilon-tilde.

r  = norm(X(1:3));
z  = X(3);
v2 = X(4:6).' * X(4:6);
V  = -mu/r + mu*J2*Re^2*(3*z^2/r^2 - 1)/(2*r^3);
E  = 0.5*v2 + V;
end

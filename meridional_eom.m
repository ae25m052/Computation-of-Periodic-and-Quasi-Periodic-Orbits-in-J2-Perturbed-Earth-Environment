function ds = meridional_eom(~, s, mu, epsilon, J2, Re, hz)
%MERIDIONAL_EOM  Reduced 2-DOF equations of motion in (rho, z).
%   ds = meridional_eom(t, s, mu, epsilon, J2, Re, hz)
%   State  s = [rho ; z ; p_rho ; p_z],  with p_rho = d(rho)/dt, p_z = dz/dt.
%       rho'' = dU/drho + hz^2 / rho^3
%       z''   = dU/dz

[~, dU] = meridional_potential(s(1), s(2), mu, epsilon, J2, Re);

ds = [ s(3) ; ...
       s(4) ; ...
       dU(1) + hz^2 / s(1)^3 ; ...
       dU(2) ];
end

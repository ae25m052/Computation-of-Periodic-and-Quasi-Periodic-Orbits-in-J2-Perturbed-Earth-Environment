function [U, dU, HU] = meridional_potential(rho, z, mu, epsilon, J2, Re)
%MERIDIONAL_POTENTIAL  Axisymmetric J2 potential in the meridional half-plane.
%   [U, dU, HU] = meridional_potential(rho, z, mu, epsilon, J2, Re)
%   U(rho,z) = mu/r - C*( 3 z^2 r^-5 - r^-3 ),   C = epsilon*J2*mu*Re^2/2
%   Returns
%     U   scalar potential
%     dU  [dU/drho ; dU/dz]                       (2x1)
%     HU  [U_rr U_rz ; U_rz U_zz]                 (2x2, exactly symmetric)

r  = sqrt(rho.^2 + z.^2);
C  = 0.5 * epsilon * J2 * mu * Re^2;

U = mu./r - C.*( 3*z.^2.*r.^-5 - r.^-3 );

if nargout > 1
    dU = [ -mu*rho.*r.^-3 + 15*C*z.^2.*rho.*r.^-7 -  3*C*rho.*r.^-5 ; ...
           -mu*z  .*r.^-3 -  9*C*z  .*r.^-5       + 15*C*z.^3.*r.^-7 ];
end

if nargout > 2
    Urr = -mu*r.^-3 + 3*mu*rho.^2.*r.^-5 ...
          + 15*C*z.^2.*r.^-7 - 105*C*z.^2.*rho.^2.*r.^-9 ...
          -  3*C*r.^-5       +  15*C*rho.^2.*r.^-7;

    Urz =  3*mu*rho.*z.*r.^-5 + 45*C*rho.*z.*r.^-7 - 105*C*rho.*z.^3.*r.^-9;

    Uzz = -mu*r.^-3 + 3*mu*z.^2.*r.^-5 ...
          -  9*C*r.^-5 + 90*C*z.^2.*r.^-7 - 105*C*z.^4.*r.^-9;

    HU = [Urr Urz ; Urz Uzz];
end
end

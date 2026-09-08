function dY = meridional_stm_eom(t, Y, mu, epsilon, J2, Re, hz)
%MERIDIONAL_STM_EOM  State + STM for the reduced system, 20 equations.
%   Y = [ s(4) ; Phi(:) (16) ]        Phi' = A(t) Phi ,  Phi(0) = I
%       A = [  0    0    1    0
%              0    0    0    1
%             a31  a32   0    0
%             a41  a42   0    0 ]
%   with
%       a31 = U_rr - 3 hz^2 / rho^4      
%       a32 = U_rz                            
%       a41 = U_rz                           
%       a42 = U_zz

s   = Y(1:4);
Phi = reshape(Y(5:20), 4, 4);

[~, dU, HU] = meridional_potential(s(1), s(2), mu, epsilon, J2, Re);

ds = [ s(3) ; s(4) ; dU(1) + hz^2/s(1)^3 ; dU(2) ];

A = zeros(4,4);
A(1,3) = 1;
A(2,4) = 1;
A(3,1) = HU(1,1) - 3*hz^2 / s(1)^4;
A(3,2) = HU(1,2);
A(4,1) = HU(2,1);
A(4,2) = HU(2,2);

dPhi = A * Phi;

dY = [ ds ; dPhi(:) ];
end

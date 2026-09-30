function dY = j2_stm_eom(t,Y,mu,Re,J2,eps)
% 42 equations: the 6 state equations plus Phidot = A*Phi.

X   = Y(1:6);
Phi = reshape(Y(7:42),6,6);

A = j2_A_matrix(X(1:3),mu,Re,J2,eps);

dY = [ j2_eom(t,X,mu,Re,J2,eps)
       reshape(A*Phi,36,1) ];

end

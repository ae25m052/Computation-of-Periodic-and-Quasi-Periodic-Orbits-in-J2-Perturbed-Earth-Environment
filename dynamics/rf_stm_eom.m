function dY = rf_stm_eom(t,Y,mu,Re,J2,Wf)
% 48 equations in the frame turning about z at rate Wf:
% state X (6), STM Phi (36), psi = dX/dWf (6).

X   = Y(1:6);
Phi = reshape(Y(7:42),6,6);
psi = Y(43:48);

% J2 part of A, then the centrifugal and Coriolis entries
A = j2_A_matrix(X(1:3),mu,Re,J2,1);
A(4,1) = A(4,1) + Wf^2;
A(5,2) = A(5,2) + Wf^2;
A(4,5) =  2*Wf;
A(5,4) = -2*Wf;

% how the equations change with Wf
dfdW = [0; 0; 0; 2*X(5) + 2*Wf*X(1); -2*X(4) + 2*Wf*X(2); 0];

dY = [ rf_eom(t,X,mu,Re,J2,Wf)
       reshape(A*Phi,36,1)
       A*psi + dfdW ];

end

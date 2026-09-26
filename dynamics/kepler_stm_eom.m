function dY = kepler_stm_eom(~, Y, mu)
%KEPLER_STM_EOM  Right-hand side for state + state transition matrix (Kepler only).
%   dY = KEPLER_STM_EOM(t, Y, mu)
%   Y is a 42x1 stacked vector:
%       Y(1:3)   position   [x; y; z]
%       Y(4:6)   velocity   [xdot; ydot; zdot]
%       Y(7:42)  the 6x6 STM Phi, unrolled column by column

r_vec = Y(1:3);
v_vec = Y(4:6);
Phi   = reshape(Y(7:42), 6, 6);    

% ---- state equations ---------------------------------------------------
r     = norm(r_vec);
a_vec = -mu * r_vec / r^3;

% ---- variational equations --------------------------------------------
A      = kepler_A_matrix(r_vec, mu);
dPhi   = A * Phi;

% ---- restack into one 42x1 column -------------------------------------
dY = [v_vec;
      a_vec;
      dPhi(:)];

end

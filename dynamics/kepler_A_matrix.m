function [A, G] = kepler_A_matrix(r_vec, mu)
%KEPLER_A_MATRIX  Jacobian A = df/dX for the unperturbed two-body problem.
%   [A, G] = KEPLER_A_MATRIX(r_vec, mu)
%   Inputs
%     r_vec : 3x1 position vector [x; y; z]
%     mu    : gravitational parameter GM (1 in canonical units)
%   Outputs
%     A : 6x6 Jacobian of the equations of motion, A = [0 I ; G 0]
%     G : 3x3 gravity-gradient block, G(i,j) = d(f_i)/d(x_j)

r_vec = r_vec(:);         

x = r_vec(1);
y = r_vec(2);
z = r_vec(3);

r  = sqrt(x*x + y*y + z*z);
r3 = r^3;
r5 = r^5;

% ---- the nine entries, written out one by one --------------------------
G = zeros(3,3);

G(1,1) = 3*mu*x*x/r5 - mu/r3;
G(1,2) = 3*mu*x*y/r5;
G(1,3) = 3*mu*x*z/r5;

G(2,1) = 3*mu*x*y/r5;              
G(2,2) = 3*mu*y*y/r5 - mu/r3;
G(2,3) = 3*mu*y*z/r5;

G(3,1) = 3*mu*x*z/r5;              
G(3,2) = 3*mu*y*z/r5;              
G(3,3) = 3*mu*z*z/r5 - mu/r3;


A = [zeros(3,3), eye(3);
     G,          zeros(3,3)];

end

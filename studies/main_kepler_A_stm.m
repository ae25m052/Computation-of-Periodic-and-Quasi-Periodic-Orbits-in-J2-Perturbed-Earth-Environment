%% MAIN_KEPLER_A_STM
%  Analytical A matrix for the unperturbed two-body problem, and the state
%  transition matrix obtained by integrating Phidot = A*Phi.
%
%  This is the epsilon = 0 case only. No J2.

clear; clc; close all;

if isempty(which('ode89'))
    error('ode89 not found. It needs MATLAB R2021b or newer.');
end

%% ---------------------------------------------------------------------
%  1. Constants and canonical units
%  ---------------------------------------------------------------------
DU = 6378.1363;                 
mu = 1;                         
mu_km = 398600.4418;            
TU = sqrt(DU^3 / mu_km);        

fprintf('Canonical units:  1 DU = %.4f km,  1 TU = %.4f s\n', DU, TU);

%% ---------------------------------------------------------------------
%  2. Seed orbit (the standard one from the project)
%  ---------------------------------------------------------------------
a_km  = 7500;
a     = a_km / DU;              
ecc   = 0.05;
inc   = deg2rad(45);
RAAN  = deg2rad(30);
argp  = deg2rad(20);
nu0   = 0;

[r0, v0] = local_coe2rv(a, ecc, inc, RAAN, argp, nu0, mu);
X0 = [r0; v0];

T = 2*pi*sqrt(a^3 / mu);        
fprintf('Seed orbit: a = %.1f km, e = %.2f, i = %.0f deg\n', a_km, ecc, rad2deg(inc));
fprintf('Period      = %.6f TU = %.5f min\n\n', T, T*TU/60);

%% ---------------------------------------------------------------------
%  3. Build A at the initial point and inspect it
%  ---------------------------------------------------------------------
[A0, G0] = kepler_A_matrix(r0, mu);

fprintf('--- A matrix at t = 0 ---\n');
disp(A0);

%% ---------------------------------------------------------------------
%  CHECK 1. G must be symmetric, and its trace must be zero
%  ---------------------------------------------------------------------
asym  = max(max(abs(G0 - G0.')));
trG   = abs(trace(G0));

fprintf('CHECK 1  symmetry  max|G - G''|  = %.3e   (expect 0)\n', asym);
fprintf('CHECK 1  Laplace   |trace(G)|    = %.3e   (expect ~1e-16)\n', trG);

%% ---------------------------------------------------------------------
%  CHECK 2. G against central finite differences of the acceleration
%  ---------------------------------------------------------------------
h    = 1e-7;
G_fd = zeros(3,3);
for j = 1:3
    ej = zeros(3,1);  ej(j) = h;
    a_plus  = local_accel(r0 + ej, mu);
    a_minus = local_accel(r0 - ej, mu);
    G_fd(:,j) = (a_plus - a_minus) / (2*h);
end
err_G = max(max(abs(G0 - G_fd)));
fprintf('CHECK 2  G vs finite differences = %.3e   (expect ~1e-9, the test''s own floor)\n\n', err_G);

%% ---------------------------------------------------------------------
%  4. Propagate state + STM over one period
%  ---------------------------------------------------------------------
tol  = 1e-13;
opts = odeset('RelTol', tol, 'AbsTol', tol);

Y0 = [X0; reshape(eye(6), 36, 1)];          
[~, Y] = ode89(@(t,Y) kepler_stm_eom(t, Y, mu), [0 T], Y0, opts);

XT  = Y(end, 1:6).';
Phi = reshape(Y(end, 7:42).', 6, 6);

fprintf('--- STM Phi(T) ---\n');
disp(Phi);

%% ---------------------------------------------------------------------
%  CHECK 3. Phi against finite differences of the flow
%  ---------------------------------------------------------------------
Phi_fd = zeros(6,6);
for j = 1:6
    ej = zeros(6,1);  ej(j) = h;
    Xp = local_prop_state(X0 + ej, T, mu, opts);
    Xm = local_prop_state(X0 - ej, T, mu, opts);
    Phi_fd(:,j) = (Xp - Xm) / (2*h);
end
err_Phi = max(max(abs(Phi - Phi_fd)));
fprintf('CHECK 3  Phi vs finite differences = %.3e   (expect ~1e-8, the test''s floor at h = 1e-7)\n', err_Phi);

%% ---------------------------------------------------------------------
%  CHECK 4. Symplecticity
%  ---------------------------------------------------------------------
J = [zeros(3,3), eye(3);
     -eye(3),    zeros(3,3)];
symp = max(max(abs(Phi.' * J * Phi - J)));
fprintf('CHECK 4  symplecticity |Phi''*J*Phi - J| = %.3e   (expect ~1e-14)\n', symp);

%% ---------------------------------------------------------------------
%  CHECK 5. Determinant
%  ---------------------------------------------------------------------
fprintf('CHECK 5  |det(Phi) - 1| = %.3e   (expect ~1e-14)\n', abs(det(Phi) - 1));

%% ---------------------------------------------------------------------
%  CHECK 6. Closure of the orbit
%  ---------------------------------------------------------------------
closure_km = norm(XT(1:3) - X0(1:3)) * DU;
fprintf('CHECK 6  closure after one period = %.3e km   (expect ~1e-9 km)\n\n', closure_km);

%% ---------------------------------------------------------------------
%  CHECK 7. Eigenvalues, and the structure of Phi - I
%  ---------------------------------------------------------------------
ev = eig(Phi);
fprintf('CHECK 7  eigenvalues of Phi(T):\n');
for k = 1:6
    fprintf('           %+.12f %+.3e i\n', real(ev(k)), imag(ev(k)));
end
fprintf('CHECK 7  max | |lambda| - 1 | = %.3e   (expect ~1e-13)\n', max(abs(abs(ev) - 1)));

sv = svd(Phi - eye(6));
fprintf('CHECK 7  singular values of Phi - I:\n           ');
fprintf('%.3e  ', sv);
fprintf('\n           one large value, five at integrator noise: rank 1, as expected.\n');

fprintf('CHECK 7  max|(Phi - I)^2| = %.3e   (expect small: shear, so the square is zero)\n', ...
        max(max(abs((Phi - eye(6))^2))));

fprintf('\nDone.\n');


%% =====================================================================
%  Local functions. In MATLAB these must sit at the very end of the script.
%  =====================================================================

function a_vec = local_accel(r_vec, mu)
r_vec = r_vec(:);
a_vec = -mu * r_vec / norm(r_vec)^3;
end


function XT = local_prop_state(X0, T, mu, opts)
f = @(t, X) [X(4:6); -mu*X(1:3)/norm(X(1:3))^3];
[~, X] = ode89(f, [0 T], X0(:), opts);
XT = X(end,:).';
end


function [r_vec, v_vec] = local_coe2rv(a, ecc, inc, RAAN, argp, nu, mu)
p = a * (1 - ecc^2);
r = p / (1 + ecc*cos(nu));

r_pf = [r*cos(nu); r*sin(nu); 0];
v_pf = sqrt(mu/p) * [-sin(nu); ecc + cos(nu); 0];

cO = cos(RAAN); sO = sin(RAAN);
cw = cos(argp); sw = sin(argp);
ci = cos(inc);  si = sin(inc);

R = [cO*cw - sO*sw*ci,  -cO*sw - sO*cw*ci,   sO*si;
     sO*cw + cO*sw*ci,  -sO*sw + cO*cw*ci,  -cO*si;
     sw*si,              cw*si,              ci   ];

r_vec = R * r_pf;
v_vec = R * v_pf;
end

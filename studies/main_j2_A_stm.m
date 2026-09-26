%% MAIN_J2_A_STM
%  Analytical A matrix for the two-body + J2 problem, and the state
%  transition matrix obtained by integrating Phidot = A*Phi.

clear; clc; close all;

if isempty(which('ode89'))
    error('ode89 not found. It needs MATLAB R2021b or newer.');
end

%% ---------------------------------------------------------------------
%  1. Constants and canonical units
%  ---------------------------------------------------------------------
DU     = 6378.1363;                
mu_km  = 398600.4418;               
TU     = sqrt(DU^3 / mu_km);       
mu     = 1;                         
Re     = 1;                         
J2     = 1.08262668e-3;

fprintf('Canonical units:  1 DU = %.4f km,  1 TU = %.4f s\n', DU, TU);

%% ---------------------------------------------------------------------
%  2. Seed orbit
%  ---------------------------------------------------------------------
a_km = 7500;
a    = a_km / DU;
ecc  = 0.05;
inc  = deg2rad(45);
RAAN = deg2rad(30);
argp = deg2rad(20);
nu0  = 0;

[r0, v0] = local_coe2rv(a, ecc, inc, RAAN, argp, nu0, mu);
X0 = [r0; v0];

T = 2*pi*sqrt(a^3 / mu);            
fprintf('Seed orbit: a = %.1f km, e = %.2f, i = %.0f deg\n', a_km, ecc, rad2deg(inc));
fprintf('Reference time T = %.6f TU = %.5f min\n', T, T*TU/60);
fprintf(['(At eps = 1 this is NOT a period. It is just a fixed span over\n' ...
         ' which to propagate, so the two cases can be compared.)\n']);

tol  = 1e-13;
opts = odeset('RelTol', tol, 'AbsTol', tol);
h    = 1e-7;                       

J = [zeros(3,3), eye(3);
     -eye(3),    zeros(3,3)];       

%% ---------------------------------------------------------------------
%  3. Run every check at eps = 0 and eps = 1
%  ---------------------------------------------------------------------
eps_list = [0, 1];

for ii = 1:numel(eps_list)

    eps = eps_list(ii);
    fprintf('\n=====================================================\n');
    fprintf('  eps = %g\n', eps);
    fprintf('=====================================================\n');

    % ---- build A at t = 0 ---------------------------------------------
    [A0, G0] = j2_A_matrix(r0, mu, eps, J2, Re);
    fprintf('\nG at t = 0:\n');
    disp(G0);

    % ---- CHECK 1: symmetry and zero trace ------------------------------
    fprintf('CHECK 1  max|G - G''|  = %.3e   (expect 0)\n', max(max(abs(G0-G0.'))));
    fprintf('CHECK 1  |trace(G)|    = %.3e   (expect ~1e-16)\n', abs(trace(G0)));

    % ---- CHECK 2: G against finite differences of the acceleration -----
    G_fd = zeros(3,3);
    for j = 1:3
        ej = zeros(3,1);  ej(j) = h;
        G_fd(:,j) = (local_accel(r0+ej, mu, eps, J2, Re) - ...
                     local_accel(r0-ej, mu, eps, J2, Re)) / (2*h);
    end
    fprintf('CHECK 2  G vs finite differences = %.3e   (expect ~1e-9)\n', ...
            max(max(abs(G0 - G_fd))));

    % ---- propagate state + STM ----------------------------------------
    Y0 = [X0; reshape(eye(6), 36, 1)];
    [~, Y] = ode89(@(t,Y) j2_stm_eom(t, Y, mu, eps, J2, Re), [0 T], Y0, opts);
    XT  = Y(end, 1:6).';
    Phi = reshape(Y(end, 7:42).', 6, 6);

    fprintf('\nPhi(T):\n');
    disp(Phi);

    % ---- CHECK 3: Phi against finite differences of the flow -----------
    Phi_fd = zeros(6,6);
    for j = 1:6
        ej = zeros(6,1);  ej(j) = h;
        Xp = local_prop_state(X0+ej, T, mu, eps, J2, Re, opts);
        Xm = local_prop_state(X0-ej, T, mu, eps, J2, Re, opts);
        Phi_fd(:,j) = (Xp - Xm) / (2*h);
    end
    fprintf('CHECK 3  Phi vs finite differences = %.3e   (expect ~1e-8)\n', ...
            max(max(abs(Phi - Phi_fd))));

    % ---- CHECK 4: symplecticity ---------------------------------------
    fprintf('CHECK 4  |Phi''*J*Phi - J| = %.3e   (expect ~1e-13)\n', ...
            max(max(abs(Phi.'*J*Phi - J))));

    % ---- CHECK 5: determinant -----------------------------------------
    fprintf('CHECK 5  |det(Phi) - 1| = %.3e   (expect ~1e-13)\n', abs(det(Phi)-1));

    % ---- CHECK 6: closure / drift -------------------------------------
    drift_km = norm(XT(1:3) - X0(1:3)) * DU;
    if eps == 0
        fprintf('CHECK 6  closure after one period = %.3e km   (expect ~1e-9 km)\n', drift_km);
    else
        fprintf('CHECK 6  J2 drift after one lap   = %.4f km   (physical, not an error)\n', drift_km);
    end

    % ---- CHECK 7: conserved quantities --------------------------------
    E0 = 0.5*(v0.'*v0) - local_potential(r0, mu, eps, J2, Re);
    ET = 0.5*(XT(4:6).'*XT(4:6)) - local_potential(XT(1:3), mu, eps, J2, Re);
    hz0 = r0(1)*v0(2) - r0(2)*v0(1);
    hzT = XT(1)*XT(5) - XT(2)*XT(4);
    fprintf('CHECK 7  energy drift = %.3e   (expect ~1e-15)\n', abs(ET-E0));
    fprintf('CHECK 7  h_z drift    = %.3e   (expect ~1e-15)\n', abs(hzT-hz0));

    % ---- CHECK 8: eigenvalues -----------------------------------------
    ev = eig(Phi);
    fprintf('CHECK 8  eigenvalues of Phi(T):\n');
    for m = 1:6
        fprintf('           %+.10f %+.3e i    |lambda| = %.10f\n', ...
                real(ev(m)), imag(ev(m)), abs(ev(m)));
    end
    if eps == 0
        fprintf('CHECK 8  max | |lambda| - 1 | = %.3e   (all six must be 1: shear)\n', ...
                max(abs(abs(ev)-1)));
    else
        fprintf('CHECK 8  they are no longer all 1. That is the physics of J2 appearing.\n');
    end

    % ---- CHECK 9: reciprocal pairing ----------------------------------
  [~, idx] = sort(abs(ev));
    ev_s = ev(idx);
    fprintf('CHECK 9  reciprocal pairs, lambda * (its partner) :\n');
    for m = 1:3
        fprintf('           %.12f\n', abs(ev_s(m)*ev_s(7-m)));
    end
    fprintf('           (all three must be 1)\n');

    % ---- singular values of Phi - I ------------------------------------
    sv = svd(Phi - eye(6));
    fprintf('         singular values of Phi - I:\n           ');
    fprintf('%.3e  ', sv);
    fprintf('\n');
    if eps == 0
        fprintf('           rank 1: five are at noise. This is why the corrector\n');
        fprintf('           cannot be started at eps = 0.\n');
    else
        fprintf('           all six are non-zero now, and they shrink linearly as\n');
        fprintf('           eps is reduced.\n');
    end

end

fprintf('\nDone.\n');


%% =====================================================================
%  Local functions
%  =====================================================================

function a_vec = local_accel(r_vec, mu, eps, J2, Re)
%LOCAL_ACCEL  Two-body + J2 acceleration.
r_vec = r_vec(:);
x = r_vec(1); y = r_vec(2); z = r_vec(3);
r = sqrt(x*x + y*y + z*z);
k = 1.5 * eps * J2 * mu * Re^2;
a_vec = [ -mu*x/r^3 - k*( x/r^5 - 5*x*z^2/r^7 );
          -mu*y/r^3 - k*( y/r^5 - 5*y*z^2/r^7 );
          -mu*z/r^3 - k*( 3*z/r^5 - 5*z^3/r^7 ) ];
end


function U = local_potential(r_vec, mu, eps, J2, Re)
%LOCAL_POTENTIAL  Potential with the sign convention a = +grad(U).
r_vec = r_vec(:);
z = r_vec(3);
r = norm(r_vec);
U = mu/r - eps*0.5*J2*mu*Re^2*( 3*z^2/r^5 - 1/r^3 );
end


function XT = local_prop_state(X0, T, mu, eps, J2, Re, opts)
%LOCAL_PROP_STATE  Propagate the 6 state equations only (no STM).
f = @(t, X) [X(4:6); local_accel(X(1:3), mu, eps, J2, Re)];
[~, X] = ode89(f, [0 T], X0(:), opts);
XT = X(end,:).';
end


function [r_vec, v_vec] = local_coe2rv(a, ecc, inc, RAAN, argp, nu, mu)
%LOCAL_COE2RV  Classical orbital elements to position and velocity.
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

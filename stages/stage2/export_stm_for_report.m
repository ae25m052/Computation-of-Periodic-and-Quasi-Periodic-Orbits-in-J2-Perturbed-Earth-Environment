function export_stm_for_report()
%EXPORT_STM_FOR_REPORT  Compute Phi(T) at eps = 0 and eps = 1 and print
%                       everything needed for the Stage 2 report section.
%
%   Just run it:
%       >> export_stm_for_report
%
%   It is DELIBERATELY SELF-CONTAINED: it does not call coe2rv, stm_eom,
%   gravity_gradient, params_earth or params_seed. That way it runs no
%   matter how your own function signatures are arranged, and it doubles
%   as an independent cross-check on your Stage 2 code. Your numbers and
%   these should agree to ~1e-13.
%
%   It writes two files into the current folder:
%       stm_report_output.txt   <- the whole printout (send me this)
%       stm_latex_blocks.tex    <- just the LaTeX matrices
%
%   Author aid for: Computation of Periodic and Quasi-Periodic Orbits in a
%   J2-Perturbed Earth Environment, R. P. Khismatrao (AE25M052), IIT Madras

clc;

%% ---------------- constants (must match params_earth.m) --------------
mu = 398600.4418;        % km^3/s^2
Re = 6378.1363;          % km
J2 = 1.08262668e-3;      % -

DU = Re;
TU = sqrt(DU^3/mu);
VU = DU/TU;

%% ---------------- seed orbit (must match params_seed.m) --------------
a    = 7500;             % km
ecc  = 0.05;
inc  = deg2rad(45);
RAAN = deg2rad(30);
argp = deg2rad(20);
nu   = deg2rad(0);

T = 2*pi*sqrt(a^3/mu);   % Keplerian period, s

x0 = elements_to_state(a, ecc, inc, RAAN, argp, nu, mu);

%% ---------------- integrate state + STM ------------------------------
opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

results = struct('eps', {0, 1}, 'Phi', {[], []}, 'Phi_can', {[], []});

for k = 1:2
    eps = results(k).eps;

    Y0 = [x0; reshape(eye(6), 36, 1)];
    [~, Y] = ode89(@(t,Y) stm_rhs(t, Y, mu, Re, J2, eps), [0 T], Y0, opts);

    Phi = reshape(Y(end, 7:42).', 6, 6);

    % convert to canonical units by the similarity transform Phi' = S\Phi*S
    S = diag([DU DU DU VU VU VU]);
    Phi_can = S \ Phi * S;

    results(k).Phi     = Phi;
    results(k).Phi_can = Phi_can;
end

%% ---------------- report ---------------------------------------------
fid = fopen('stm_report_output.txt', 'w');
fidtex = fopen('stm_latex_blocks.tex', 'w');
out = @(varargin) both(fid, varargin{:});

out('================================================================\n');
out(' STATE TRANSITION MATRIX OVER ONE PERIOD OF THE SEED ORBIT\n');
out('================================================================\n');
out('Seed: a = %.4f km, e = %.4f, i = %.4f deg,\n', a, ecc, rad2deg(inc));
out('      RAAN = %.4f deg, argp = %.4f deg, nu = %.4f deg\n', ...
    rad2deg(RAAN), rad2deg(argp), rad2deg(nu));
out('Period T = %.6f s = %.6f min = %.9f TU\n', T, T/60, T/TU);
out('Canonical units: DU = %.4f km, TU = %.6f s, VU = %.6f km/s\n', DU, TU, VU);
out('Integrator: ode89, RelTol = AbsTol = 1e-13\n');
out('Initial state (km, km/s):\n');
out('  r0 = [%.9f  %.9f  %.9f]\n', x0(1), x0(2), x0(3));
out('  v0 = [%.9f  %.9f  %.9f]\n', x0(4), x0(5), x0(6));

Jsym = [zeros(3) eye(3); -eye(3) zeros(3)];

for k = 1:2
    eps     = results(k).eps;
    Phi     = results(k).Phi;
    Phi_can = results(k).Phi_can;

    out('\n----------------------------------------------------------------\n');
    out(' eps = %g\n', eps);
    out('----------------------------------------------------------------\n');

    out('\nPhi(T) in PHYSICAL units (km, km/s):\n');
    print_matrix(out, Phi);

    out('\nPhi(T) in CANONICAL units (DU, VU, TU):\n');
    print_matrix(out, Phi_can);

    out('\nChecks (canonical):\n');
    out('  ||Phi^T J Phi - J||   = %.6e\n', norm(Phi_can.'*Jsym*Phi_can - Jsym));
    out('  |det Phi - 1|         = %.6e\n', abs(det(Phi_can) - 1));
    out('  rank(Phi - I)         = %d\n', rank(Phi_can - eye(6), 1e-8));
    out('  ||(Phi - I)^2||       = %.6e\n', norm((Phi_can - eye(6))^2));

    ev = eig(Phi_can);
    out('\nEigenvalues of Phi(T) (canonical):\n');
    for i = 1:6
        out('  %2d:  %+.9f %+.9fi   |lambda| = %.9f\n', ...
            i, real(ev(i)), imag(ev(i)), abs(ev(i)));
    end

    % LaTeX bmatrix of the canonical matrix
    fprintf(fidtex, '%% ---- eps = %g, canonical units ----\n', eps);
    fprintf(fidtex, '%s\n\n', latex_bmatrix(Phi_can, 6));
end

out('\n================================================================\n');
out(' Files written: stm_report_output.txt, stm_latex_blocks.tex\n');
out('================================================================\n');

fclose(fid);
fclose(fidtex);

end

%% ======================================================================
function both(fid, varargin)
fprintf(varargin{:});
fprintf(fid, varargin{:});
end

%% ======================================================================
function print_matrix(out, M)
for i = 1:6
    out('  ');
    for j = 1:6
        out('%16.9g', M(i,j));
    end
    out('\n');
end
end

%% ======================================================================
function s = latex_bmatrix(M, sig)
fmt = sprintf('%%.%dg', sig);
lines = {'\begin{bmatrix}'};
for i = 1:6
    e = cell(1,6);
    for j = 1:6
        e{j} = sprintf(fmt, M(i,j));
    end
    if i < 6, tail = ' \\'; else, tail = ''; end
    lines{end+1} = ['  ' strjoin(e, ' & ') tail]; %#ok<AGROW>
end
lines{end+1} = '\end{bmatrix}';
s = strjoin(lines, newline);
end

%% ======================================================================
function dY = stm_rhs(~, Y, mu, Re, J2, eps)
%STM_RHS  42 equations: 6 for the state, 36 for the STM.
r   = Y(1:3);
v   = Y(4:6);
Phi = reshape(Y(7:42), 6, 6);

a = accel(r, mu, Re, J2, eps);
G = grav_gradient(r, mu, Re, J2, eps);

A = [zeros(3), eye(3); G, zeros(3)];

dY = [v; a; reshape(A*Phi, 36, 1)];
end

%% ======================================================================
function a = accel(rv, mu, Re, J2, eps)
x = rv(1); y = rv(2); z = rv(3);
r = norm(rv);

a = -mu*rv/r^3;

k = 1.5*J2*mu*Re^2;
f = 1/r^5 - 5*z^2/r^7;
g = 3/r^5 - 5*z^2/r^7;

a = a + eps*[-k*x*f; -k*y*f; -k*z*g];
end

%% ======================================================================
function G = grav_gradient(rv, mu, Re, J2, eps)
%GRAV_GRADIENT  Analytical d(rddot)/dr for the two-body + eps*J2 field.
%
%   G is the Hessian of a scalar potential, so it MUST be symmetric and
%   trace-free. Each off-diagonal entry is therefore computed ONCE and
%   placed in both slots. Computing G(1,2) and G(2,1) by their two
%   separate (algebraically identical) expressions would leave them a
%   few ulp apart, because x*(A*y) and y*(A*x) round differently, and
%   exact symmetry would be lost for no reason.

x = rv(1); y = rv(2); z = rv(3);
r = norm(rv);

% ---------------- two-body part: -mu/r^3 I + 3 mu/r^5 (r r^T) --------
d = -mu/r^3;
c =  3*mu/r^5;

Txx = d + c*x*x;   Tyy = d + c*y*y;   Tzz = d + c*z*z;
Txy =     c*x*y;   Txz =     c*x*z;   Tyz =     c*y*z;

% ---------------- J2 part --------------------------------------------
% a_x = -k x f,  a_y = -k y f,  a_z = -k z g
k = 1.5*J2*mu*Re^2;

f = 1/r^5 - 5*z^2/r^7;
g = 3/r^5 - 5*z^2/r^7;

fx = -5*x/r^7  + 35*z^2*x/r^9;
fy = -5*y/r^7  + 35*z^2*y/r^9;
fz = -15*z/r^7 + 35*z^3/r^9;
gz = -25*z/r^7 + 35*z^3/r^9;

Jxx = -k*(f + x*fx);   Jyy = -k*(f + y*fy);   Jzz = -k*(g + z*gz);
Jxy = -k*x*fy;         Jxz = -k*x*fz;         Jyz = -k*y*fz;
% note: Jxz could equally be written -k*z*gx, and Jyz as -k*z*gy;
% the two forms are algebraically identical, so only one is evaluated.

% ---------------- assemble symmetrically ------------------------------
% Each entry is evaluated once, then placed in both slots.
Gxx = Txx + eps*Jxx;   Gyy = Tyy + eps*Jyy;   Gzz = Tzz + eps*Jzz;
Gxy = Txy + eps*Jxy;   Gxz = Txz + eps*Jxz;   Gyz = Tyz + eps*Jyz;

G = [ Gxx, Gxy, Gxz ;
      Gxy, Gyy, Gyz ;
      Gxz, Gyz, Gzz ];

% Symmetric by construction, so this is exact, not approximate.
assert(isequal(G, G.'), 'grav_gradient: not exactly symmetric.');
end

%% ======================================================================
function x = elements_to_state(a, e, i, RAAN, argp, nu, mu)
%ELEMENTS_TO_STATE  Classical elements (angles in RADIANS) -> [r; v].
p = a*(1 - e^2);
rmag = p/(1 + e*cos(nu));

r_pqw = [rmag*cos(nu); rmag*sin(nu); 0];
v_pqw = sqrt(mu/p)*[-sin(nu); e + cos(nu); 0];

R3_W = [ cos(-RAAN) sin(-RAAN) 0; -sin(-RAAN) cos(-RAAN) 0; 0 0 1];
R1_i = [1 0 0; 0 cos(-i) sin(-i); 0 -sin(-i) cos(-i)];
R3_w = [ cos(-argp) sin(-argp) 0; -sin(-argp) cos(-argp) 0; 0 0 1];

Q = R3_W*R1_i*R3_w;

x = [Q*r_pqw; Q*v_pqw];
end
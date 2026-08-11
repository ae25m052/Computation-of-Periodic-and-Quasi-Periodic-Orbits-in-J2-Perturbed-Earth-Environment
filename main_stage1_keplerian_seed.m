clear; clc; close all;

%% ---------------------------------------------------------------
%  1. Physical constants and non-dimensionalization
%  -----------------------------------------------------------------
mu_dim = 398600.4418;      % Earth GM [km^3/s^2]
Re_dim = 6378.1363;        % Earth equatorial radius [km]
J2     = 1.08262668e-3;    % Earth's J2 (unused this stage, kept for Stage 4)

% Canonical units (DU = Re, TU chosen so that mu_nondim = 1)
DU = Re_dim;                       % length unit [km]
TU = sqrt(DU^3 / mu_dim);          % time unit [s]
VU = DU / TU;                      % velocity unit [km/s]

mu_nd = 1;                         % non-dimensional mu, by construction
Re_nd = Re_dim / DU;               % = 1, kept explicit for clarity

fprintf('--- Non-dimensionalization ---\n');
fprintf('DU = %.6f km\n', DU);
fprintf('TU = %.6f s (%.4f min)\n', TU, TU/60);
fprintf('VU = %.6f km/s\n\n', VU);

%% ---------------------------------------------------------------
%  2. Seed orbit: classical orbital elements (EDIT AS NEEDED)
%  -----------------------------------------------------------------
a_dim   = 7500;             % semi-major axis [km]
e0      = 0.05;             % eccentricity [-]
i0      = deg2rad(45);      % inclination
RAAN0   = deg2rad(30);      % RAAN
argp0   = deg2rad(20);      % argument of periapsis
nu0     = deg2rad(0);       % true anomaly at epoch

[r0_dim, v0_dim] = coe2rv(a_dim, e0, i0, RAAN0, argp0, nu0, mu_dim);

% Non-dimensionalize the initial condition
r0_nd = r0_dim / DU;
v0_nd = v0_dim / VU;
x0_nd = [r0_nd; v0_nd];

fprintf('--- Seed orbit (dimensional) ---\n');
fprintf('a = %.3f km, e = %.4f, i = %.2f deg\n', a_dim, e0, rad2deg(i0));
fprintf('r0 = [%.6f %.6f %.6f] km\n', r0_dim);
fprintf('v0 = [%.6f %.6f %.6f] km/s\n\n', v0_dim);

%% ---------------------------------------------------------------
%  3. Propagate over several orbital periods
%  -----------------------------------------------------------------
T_dim = 2*pi*sqrt(a_dim^3 / mu_dim);   % analytic Keplerian period [s]
T_nd  = T_dim / TU;                     % non-dimensional period

n_periods = 3;
tspan = linspace(0, n_periods*T_nd, 2000);

opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

% epsilon = 0 => pure Keplerian two-body.
% Set >0 later for J2 continuation.
epsilon = 0;

[t_nd, X_nd] = ode89(@(t,x) twobody_eom(t, x, mu_nd, epsilon, J2, Re_nd), ...
                     tspan, x0_nd, opts);
X_nd = X_nd.';   % 6xN, columns = states, to match compute_invariants

fprintf('--- Propagation ---\n');
fprintf('Integrator: ode89, RelTol = AbsTol = 1e-13\n');
fprintf('Analytic Keplerian period T = %.4f s (%.4f TU)\n', T_dim, T_nd);
fprintf('Propagated for %d periods, %d output samples\n\n', n_periods, numel(t_nd));

%% ---------------------------------------------------------------
%  4. Validation Check 1-3: Energy, Angular Momentum, LRL vector
%  -----------------------------------------------------------------
[energy_nd, h_vec_nd, e_vec_nd] = compute_invariants(X_nd, mu_nd);

energy0 = energy_nd(1);
h0_vec  = h_vec_nd(:,1);
e0_vec  = e_vec_nd(:,1);

energy_drift = abs(energy_nd - energy0);
h_drift      = vecnorm(h_vec_nd - h0_vec);
e_drift      = vecnorm(e_vec_nd - e0_vec);

fprintf('--- Conservation checks (non-dimensional) ---\n');
fprintf('Specific energy at t0        : %.12e\n', energy0);
fprintf('Max |energy(t) - energy(0)|  : %.3e\n', max(energy_drift));
fprintf('Ang. momentum |h| at t0      : %.12e\n', norm(h0_vec));
fprintf('Max |h(t) - h(0)| (vector)   : %.3e\n', max(h_drift));
fprintf('Eccentricity |e| at t0       : %.12e  (expect %.6f)\n', norm(e0_vec), e0);
fprintf('Max |e_vec(t) - e_vec(0)|    : %.3e\n\n', max(e_drift));

%% ---------------------------------------------------------------
%  5. Validation Check 4: Closure error (state after one period)
%  -----------------------------------------------------------------
% Re-integrate exactly one period and compare final state to initial.
[~, X_close] = ode89(@(t,x) twobody_eom(t, x, mu_nd, epsilon, J2, Re_nd), ...
                      [0 T_nd], x0_nd, opts);
xf_nd = X_close(end, :).';

closure_error_r = norm(xf_nd(1:3) - x0_nd(1:3));
closure_error_v = norm(xf_nd(4:6) - x0_nd(4:6));

fprintf('--- Closure error over one period ---\n');
fprintf('|r(T) - r(0)| = %.3e (non-dim)  = %.3e km\n', closure_error_r, closure_error_r*DU);
fprintf('|v(T) - v(0)| = %.3e (non-dim)  = %.3e km/s\n\n', closure_error_v, closure_error_v*VU);

%% ---------------------------------------------------------------
%  6. Plots (for the "Initial Results" slide)
%  -----------------------------------------------------------------
% 6a. 3D orbit
orbitColor = [0.95 0.75 0.05];   % single trace colour -- edit to taste

figure('Name', 'Stage 1: Keplerian Seed Orbit', 'Color', 'w');
plot_earth(Re_nd); hold on;
plot3(X_nd(1,:), X_nd(2,:), X_nd(3,:), '-', ...
      'Color', orbitColor, 'LineWidth', 1.8);
h_init = plot3(x0_nd(1), x0_nd(2), x0_nd(3), 'w*', 'MarkerSize', 12, 'LineWidth', 1.5);
grid on; set(gca, 'GridAlpha', 0.25);
xlabel('x [DU]'); ylabel('y [DU]'); zlabel('z [DU]');
title('Two-Body Keplerian Seed Orbit (ECI, non-dim)');
legend(h_init, 'Initial state', 'Location', 'best');
view(45, 25); axis vis3d;

lim = 1.15 * max(vecnorm(X_nd(1:3,:)));
xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);
style_light(gcf, 12);

% 6b. Conservation drift vs time
figure('Name', 'Stage 1: Conservation Diagnostics', 'Color', 'w');
subplot(3,1,1);
semilogy(t_nd, energy_drift + eps, 'LineWidth', 1.2);
ylabel('|\Delta E|'); grid on; title('Specific Energy Drift');

subplot(3,1,2);
semilogy(t_nd, h_drift + eps, 'LineWidth', 1.2);
ylabel('|\Delta h|'); grid on; title('Angular Momentum Vector Drift');

subplot(3,1,3);
semilogy(t_nd, e_drift + eps, 'LineWidth', 1.2);
ylabel('|\Delta e_{vec}|'); xlabel('time [TU]'); grid on;
title('Laplace-Runge-Lenz (Eccentricity) Vector Drift');

sgt = sgtitle('Stage 1 Validation: Conservation Quantities vs Time');
set(sgt, 'Color', [0.10 0.10 0.10], 'FontWeight', 'bold', 'FontSize', 15);
style_light(gcf, 12);
set(sgt, 'Color', [0.10 0.10 0.10]);  

%% ---------------------------------------------------------------
%  7. Summary table (copy into advisor slide / report)
%  -----------------------------------------------------------------
fprintf('=========================================================\n');
fprintf(' STAGE 1 SUMMARY (for Initial Results slide)\n');
fprintf('=========================================================\n');
fprintf(' Integrator            : ode89, RelTol=AbsTol=1e-13\n');
fprintf(' Seed orbit             : a=%.1f km, e=%.3f, i=%.1f deg\n', a_dim, e0, rad2deg(i0));
fprintf(' Periods propagated    : %d\n', n_periods);
fprintf(' Max energy drift      : %.3e\n', max(energy_drift));
fprintf(' Max ang. mom. drift   : %.3e\n', max(h_drift));
fprintf(' Max LRL vector drift  : %.3e\n', max(e_drift));
fprintf(' Closure error (r)     : %.3e km\n', closure_error_r*DU);
fprintf(' Closure error (v)     : %.3e km/s\n', closure_error_v*VU);
fprintf('=========================================================\n');

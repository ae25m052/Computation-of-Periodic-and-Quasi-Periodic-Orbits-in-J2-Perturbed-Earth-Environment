clear; clc; close all;

%% ---------------------------------------------------------------
%  1. Constants and non-dimensionalization (same convention as Stage 1)
%  -----------------------------------------------------------------
mu_dim = 398600.4418;       % Earth GM [km^3/s^2]
Re_dim = 6378.1363;         % Earth equatorial radius [km]
J2     = 1.08262668e-3;     % Earth's J2 (dimensionless)

DU = Re_dim;
TU = sqrt(DU^3/mu_dim);
VU = DU/TU;

mu_nd = 1;
Re_nd = Re_dim/DU;          % = 1

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

%% ---------------------------------------------------------------
%  2. Initial condition (edit as needed)
%  -----------------------------------------------------------------
a_dim   = 7500;              % km
e0      = 0.05;
i0      = deg2rad(45);
RAAN0   = deg2rad(30);
argp0   = deg2rad(20);
nu0     = 0;

[r0_dim, v0_dim] = coe2rv(a_dim, e0, i0, RAAN0, argp0, nu0, mu_dim);
r0_nd = r0_dim/DU;  v0_nd = v0_dim/VU;
x0_nd = [r0_nd; v0_nd];

T_dim = 2*pi*sqrt(a_dim^3/mu_dim);
T_nd  = T_dim/TU;

fprintf('Seed orbit: a=%.1f km, e=%.3f, i=%.1f deg, period T=%.2f min\n\n', ...
        a_dim, e0, rad2deg(i0), T_dim/60);

%% ---------------------------------------------------------------
%  PART A: Pure Keplerian propagation (epsilon = 0) -> CLOSED orbit
%  -----------------------------------------------------------------
n_periods_A = 3;
tspan_A = linspace(0, n_periods_A*T_nd, 1500);

[t_A, X_A] = ode89(@(t,x) twobody_eom(t, x, mu_nd, 0, J2, Re_nd), ...
                    tspan_A, x0_nd, opts);
X_A = X_A.';

% Closure check: state after exactly one period
[~, Xc_A] = ode89(@(t,x) twobody_eom(t, x, mu_nd, 0, J2, Re_nd), ...
                   [0 T_nd], x0_nd, opts);
xf_A = Xc_A(end,:).';
closure_r_A = norm(xf_A(1:3) - x0_nd(1:3)) * DU;   % km

fprintf('--- PART A: Keplerian (epsilon = 0) ---\n');
fprintf('Closure error after 1 period: %.3e km  (orbit closes)\n\n', closure_r_A);

%% ---------------------------------------------------------------
%  PART B: J2-perturbed propagation (epsilon = 1) -> orbit precesses
%  -----------------------------------------------------------------
n_periods_B = 30;    % propagate many periods so precession is visible
tspan_B = linspace(0, n_periods_B*T_nd, 6000);

[t_B, X_B] = ode89(@(t,x) twobody_eom(t, x, mu_nd, 1, J2, Re_nd), ...
                    tspan_B, x0_nd, opts);
X_B = X_B.';

[~, Xc_B] = ode89(@(t,x) twobody_eom(t, x, mu_nd, 1, J2, Re_nd), ...
                   [0 T_nd], x0_nd, opts);
xf_B = Xc_B(end,:).';
closure_r_B = norm(xf_B(1:3) - x0_nd(1:3)) * DU;   % km

fprintf('--- PART B: J2-perturbed (epsilon = 1) ---\n');
fprintf('"Closure" error after 1 nominal Keplerian period: %.3e km\n', closure_r_B);
fprintf('(large and non-decaying with more periods => orbit is NOT closed)\n\n');

%% ---------------------------------------------------------------
%  3. Extract osculating elements over time, quantify precession
%  -----------------------------------------------------------------
N_B = size(X_B, 2);
RAAN_t = zeros(1, N_B);
argp_t = zeros(1, N_B);

for k = 1:N_B
    [~, ~, ~, RAAN_t(k), argp_t(k), ~] = rv2coe(X_B(1:3,k), X_B(4:6,k), mu_nd);
end

RAAN_unwrapped = unwrap(RAAN_t);
argp_unwrapped = unwrap(argp_t);

n_mean = sqrt(mu_nd/a_dim^3 * DU^3);  
a_nd = a_dim/DU;
n_mean_nd = sqrt(mu_nd/a_nd^3);               % mean motion, non-dim [1/TU]
p_nd = a_nd*(1-e0^2);

RAAN_dot_analytical = -1.5 * n_mean_nd * J2 * (Re_nd/p_nd)^2 * cos(i0);
argp_dot_analytical =  0.75 * n_mean_nd * J2 * (Re_nd/p_nd)^2 * (5*cos(i0)^2 - 1);

pfit_RAAN = polyfit(t_B, RAAN_unwrapped, 1);
pfit_argp = polyfit(t_B, argp_unwrapped, 1);

fprintf('--- J2 secular drift validation ---\n');
fprintf('RAAN_dot: analytical = %.6e /TU,  fit from propagation = %.6e /TU\n', ...
        RAAN_dot_analytical, pfit_RAAN(1));
fprintf('argp_dot: analytical = %.6e /TU,  fit from propagation = %.6e /TU\n', ...
        argp_dot_analytical, pfit_argp(1));
fprintf('Relative error (RAAN_dot): %.3e\n', ...
        abs(pfit_RAAN(1)-RAAN_dot_analytical)/abs(RAAN_dot_analytical));
fprintf('Relative error (argp_dot): %.3e\n\n', ...
        abs(pfit_argp(1)-argp_dot_analytical)/abs(argp_dot_analytical));

%% ---------------------------------------------------------------
%  4. Visualization
%  -----------------------------------------------------------------
% 4a. Side-by-side 3D orbits: closed Keplerian vs precessing J2 orbit
colorKep = [0.20 0.85 0.95];   % Part A trace colour (cyan)
colorJ2  = [1.00 0.30 0.20];   % Part B trace colour (red-orange)

figure('Name', 'Keplerian vs J2-Perturbed Orbit', 'Color', 'w', ...
       'Position', [100 100 1200 560]);

subplot(1,2,1);
plot_earth(Re_nd); hold on;
plot3(X_A(1,:), X_A(2,:), X_A(3,:), '-', 'Color', colorKep, 'LineWidth', 1.8);
grid on; set(gca, 'GridAlpha', 0.25); view(45,25); axis vis3d;
xlabel('x [DU]'); ylabel('y [DU]'); zlabel('z [DU]');
title(sprintf('Part A: Keplerian (\\epsilon=0)\n%d periods -- closed orbit', n_periods_A));
limA = 1.15 * max(vecnorm(X_A(1:3,:)));
xlim([-limA limA]); ylim([-limA limA]); zlim([-limA limA]);

subplot(1,2,2);
plot_earth(Re_nd); hold on;
plot3(X_B(1,:), X_B(2,:), X_B(3,:), '-', 'Color', colorJ2, 'LineWidth', 0.9);
grid on; set(gca, 'GridAlpha', 0.25); view(45,25); axis vis3d;
xlabel('x [DU]'); ylabel('y [DU]'); zlabel('z [DU]');
title(sprintf('Part B: J2-perturbed (\\epsilon=1)\n%d periods -- precessing (not closed)', n_periods_B));
limB = 1.15 * max(vecnorm(X_B(1:3,:)));
xlim([-limB limB]); ylim([-limB limB]); zlim([-limB limB]);

sgt = sgtitle('Governing Equations: Two-Body vs. J_2-Perturbed Motion', ...
              'FontWeight', 'bold');
set(sgt, 'Color', [0.10 0.10 0.10], 'FontSize', 15);
style_light(gcf, 12);
set(sgt, 'Color', [0.10 0.10 0.10]);   % re-assert after styling

% 4b. RAAN and argument of perigee drift vs time (numerical vs analytical)
figure('Name', 'J2 Secular Drift', 'Color', 'w');
subplot(2,1,1);
plot(t_B, rad2deg(RAAN_unwrapped), 'Color', [0.85 0.2 0.2], 'LineWidth', 1.4); hold on;
plot(t_B, rad2deg(polyval(pfit_RAAN, t_B)), 'k--', 'LineWidth', 1);
ylabel('RAAN [deg]'); grid on; set(gca, 'GridAlpha', 0.3);
legend('Propagated (osculating)', 'Linear fit', 'Location', 'best');
title('Nodal Regression: d\Omega/dt due to J_2');

subplot(2,1,2);
plot(t_B, rad2deg(argp_unwrapped), 'Color', [0.85 0.2 0.2], 'LineWidth', 1.4); hold on;
plot(t_B, rad2deg(polyval(pfit_argp, t_B)), 'k--', 'LineWidth', 1);
xlabel('time [TU]'); ylabel('\omega [deg]'); grid on; set(gca, 'GridAlpha', 0.3);
legend('Propagated (osculating)', 'Linear fit', 'Location', 'best');
title('Apsidal Drift: d\omega/dt due to J_2');
style_light(gcf, 12);

% 4c. Optional animated GIFs for the slide deck (comment out if not needed)
animate_orbit(X_A, Re_nd, 'Keplerian Orbit (no J2)', 'keplerian_orbit.gif', 3);
animate_orbit(X_B, Re_nd, 'J2-Perturbed Orbit (precessing)', 'j2_orbit.gif', 8);

fprintf('=========================================================\n');
fprintf(' SUMMARY (for slide: Keplerian vs J2)\n');
fprintf('=========================================================\n');
fprintf(' Part A closure error (no J2)  : %.3e km  -> orbit closes\n', closure_r_A);
fprintf(' Part B "closure" error (J2 on): %.3e km  -> orbit does NOT close\n', closure_r_B);
fprintf(' RAAN drift rate  (num vs analytical): %.3e vs %.3e /TU\n', pfit_RAAN(1), RAAN_dot_analytical);
fprintf(' argp drift rate  (num vs analytical): %.3e vs %.3e /TU\n', pfit_argp(1), argp_dot_analytical);
fprintf('=========================================================\n');

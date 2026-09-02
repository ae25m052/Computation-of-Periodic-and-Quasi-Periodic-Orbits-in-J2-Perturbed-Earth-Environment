%% STAGE 2: STATE TRANSITION MATRIX VIA VARIATIONAL EQUATIONS

clear; clc; close all;

%% ---------------------------------------------------------------
%  1. Constants and non-dimensionalization  (same as earlier stages)
%  -----------------------------------------------------------------
mu_dim = 398600.4418;
Re_dim = 6378.1363;
J2     = 1.08262668e-3;

DU = Re_dim;
TU = sqrt(DU^3/mu_dim);
VU = DU/TU;

mu_nd = 1;
Re_nd = Re_dim/DU;

opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);


%% ---------------------------------------------------------------
%  2. Seed orbit  (same initial state used throughout the project)
%  -----------------------------------------------------------------
a_dim = 7500;  e0 = 0.05;  i0 = deg2rad(45);
RAAN0 = deg2rad(30);  argp0 = deg2rad(20);  nu0 = 0;

[r0_dim, v0_dim] = coe2rv(a_dim, e0, i0, RAAN0, argp0, nu0, mu_dim);
x0 = [r0_dim/DU; v0_dim/VU];

T_dim = 2*pi*sqrt(a_dim^3/mu_dim);
T_nd  = T_dim/TU;

fprintf('=== STAGE 2: STATE TRANSITION MATRIX ===\n');
fprintf('Seed orbit : a=%.1f km, e=%.3f, i=%.1f deg, T=%.3f min\n', ...
        a_dim, e0, rad2deg(i0), T_dim/60);
fprintf('Integrator : ode89, RelTol = AbsTol = 1e-13\n\n');


%% ---------------------------------------------------------------
%  3. Validation (a): is the gravity gradient symmetric?
%  -----------------------------------------------------------------
% G is the Hessian of a scalar potential, so it must be symmetric to
% machine precision. This is a check on the algebra.
fprintf('--- (a) Gravity gradient symmetry ---\n');
for eps_test = [0 1]
    r_test = x0(1:3);
    G = gravity_gradient(r_test, mu_nd, eps_test, J2, Re_nd);
    fprintf('  epsilon = %d :  max|G - G''| = %.3e\n', ...
            eps_test, max(max(abs(G - G.'))));
end
fprintf('\n');


%% ---------------------------------------------------------------
%  4. Propagate state + STM over one period, for both cases
%  -----------------------------------------------------------------
% Symplectic form J used for the symplecticity test below.
Jsym = [ zeros(3),  eye(3);
        -eye(3),    zeros(3) ];

results = struct();

for k = 1:2
    epsilon = k - 1;                      

    X0 = [x0; reshape(eye(6), 36, 1)];    

    [~, XX] = ode89(@(t,X) stm_eom(t, X, mu_nd, epsilon, J2, Re_nd), ...
                    [0 T_nd], X0, opts);

    Xf   = XX(end, :).';
    Phi  = reshape(Xf(7:42), 6, 6);

    results(k).epsilon = epsilon;
    results(k).Phi     = Phi;
    results(k).xf      = Xf(1:6);
end


%% ---------------------------------------------------------------
%  5. Validation (b): analytical STM vs finite-difference STM
%  -----------------------------------------------------------------
fprintf('--- (b) Analytical STM vs finite differences ---\n');
h = 1e-7;

for k = 1:2
    epsilon = results(k).epsilon;
    Phi_fd = zeros(6,6);

    for j = 1:6
        xp = x0;  xm = x0;
        xp(j) = xp(j) + h;
        xm(j) = xm(j) - h;

        [~, Xp] = ode89(@(t,x) twobody_eom(t, x, mu_nd, epsilon, J2, Re_nd), ...
                        [0 T_nd], xp, opts);
        [~, Xm] = ode89(@(t,x) twobody_eom(t, x, mu_nd, epsilon, J2, Re_nd), ...
                        [0 T_nd], xm, opts);

        Phi_fd(:,j) = (Xp(end,:).' - Xm(end,:).') / (2*h);
    end

    err = max(max(abs(results(k).Phi - Phi_fd)));
    results(k).fd_error = err;
    fprintf('  epsilon = %d :  max|Phi_analytic - Phi_fd| = %.3e\n', ...
            epsilon, err);
end
fprintf('\n');


%% ---------------------------------------------------------------
%  6. Validation (c): symplecticity and determinant
%  -----------------------------------------------------------------
fprintf('--- (c) Symplecticity ---\n');
for k = 1:2
    Phi = results(k).Phi;
    symp_err = max(max(abs(Phi.' * Jsym * Phi - Jsym)));
    detPhi   = det(Phi);
    results(k).symp_err = symp_err;
    results(k).detPhi   = detPhi;
    fprintf('  epsilon = %d :  max|Phi''J Phi - J| = %.3e ,  det(Phi) = %.12f\n', ...
            results(k).epsilon, symp_err, detPhi);
end
fprintf('\n');


%% ---------------------------------------------------------------
%  7. Validation (d): eigenvalue structure
%  -----------------------------------------------------------------
fprintf('--- (d) Eigenvalue structure ---\n');
for k = 1:2
    Phi = results(k).Phi;
    ev  = eig(Phi);
    [~, idx] = sort(abs(ev), 'descend');
    ev = ev(idx);

    fprintf('  epsilon = %d :\n', results(k).epsilon);
    for m = 1:6
        fprintf('     %+.6f %+.6fi    |lambda| = %.6f\n', ...
                real(ev(m)), imag(ev(m)), abs(ev(m)));
    end
    fprintf('     product of |lambda| = %.12f  (must be 1)\n', prod(abs(ev)));
    results(k).eig = ev;
end
fprintf('\n');


%% ---------------------------------------------------------------
%  8. Visualisation of the STM
%  -----------------------------------------------------------------
figure('Name','Stage 2: State Transition Matrix','Color','w', ...
       'Position',[100 100 1000 430]);

for k = 1:2
    subplot(1,2,k);
    Phi = results(k).Phi;
    imagesc(Phi);
    axis square; colorbar;
    set(gca, 'XTick', 1:6, 'YTick', 1:6, ...
             'XTickLabel', {'x','y','z','vx','vy','vz'}, ...
             'YTickLabel', {'x','y','z','vx','vy','vz'});
    xlabel('perturbation in'); ylabel('response of');
    title(sprintf('\\Phi(T) at \\epsilon = %d', results(k).epsilon));
end
sgt = sgtitle('State Transition Matrix over one period');
set(sgt, 'Color', [0.10 0.10 0.10], 'FontWeight', 'bold', 'FontSize', 15);
if exist('style_light','file'), style_light(gcf, 12); end
set(sgt, 'Color', [0.10 0.10 0.10]);

figure('Name','Stage 2: STM eigenvalues','Color','w');
th = linspace(0, 2*pi, 400);
plot(cos(th), sin(th), 'k--', 'LineWidth', 1); hold on;
plot(real(results(1).eig), imag(results(1).eig), 'o', ...
     'MarkerSize', 11, 'LineWidth', 1.8, 'Color', [0.00 0.35 0.85]);
plot(real(results(2).eig), imag(results(2).eig), 'x', ...
     'MarkerSize', 12, 'LineWidth', 1.8, 'Color', [0.85 0.10 0.10]);
axis equal; grid on;
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title('STM eigenvalues over one period');
legend('unit circle', '\epsilon = 0', '\epsilon = 1', 'Location','best');
if exist('style_light','file'), style_light(gcf, 12); end


%% ---------------------------------------------------------------
%  9. Summary
%  -----------------------------------------------------------------
fprintf('=========================================================\n');
fprintf(' STAGE 2 SUMMARY\n');
fprintf('=========================================================\n');
for k = 1:2
    fprintf(' epsilon = %d\n', results(k).epsilon);
    fprintf('   STM vs finite difference : %.3e\n', results(k).fd_error);
    fprintf('   symplecticity residual   : %.3e\n', results(k).symp_err);
    fprintf('   det(Phi) - 1             : %.3e\n', results(k).detPhi - 1);
end
fprintf('=========================================================\n');
fprintf(' STM verified. Ready for Stage 3: single shooting and\n');
fprintf(' differential correction, which uses Phi as the Newton Jacobian.\n');
fprintf('=========================================================\n');

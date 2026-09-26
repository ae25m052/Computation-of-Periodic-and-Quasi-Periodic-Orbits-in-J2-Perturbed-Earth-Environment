%% main_dc_equatorial.m
%  Differential correction for the PROJECT orbit in the equatorial plane,
%  under the full J2 model, in ECI.
%     perigee altitude 300 km,  e = 0.65,  i = 0,  RAAN = 0,  argp = 0,  nu = 0
%  CANONICAL UNITS THROUGHOUT: DU = Re, TU = sqrt(Re^3/mu), mu = 1, Re = 1.
%  The only dimensional input (300 km) is normalised inside
%  dc_seed_perigee.m before anything else uses it. km, s and m/s appear
%  only in printed output and on some plot axes, for reading.
%
%  What it does
%    1. Kepler starting state from the elements
%    2. Differential correction: fix perigee, correct perigee speed so the
%       J2 orbit reaches the intended apogee (e_J2 = 0.65 geometrically)
%    3. One radial period: closure checks and the monodromy matrix in the
%       frame that turns with the ellipse
%    4. 400 radial periods: apsidal drift, energy and h conservation
%    5. Figures (5 PNGs)
%    6. Save CSV + MAT

clear; clc; close all;

%% ---------------------------------------------------------------------
%  0. Settings
%  ---------------------------------------------------------------------
mu = 1;  Re = 1;                   % canonical
J2 = 1.08262668e-3;                % project value (j2_constants.m)

alt_p_km = 300;                    % perigee altitude
e_orb    = 0.65;                   % eccentricity
% i = 0, RAAN = 0, argp = 0, nu = 0  (start at perigee on the +x axis)

N_REVS = 400;                      % long run, as in the inclination sweep
opts   = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

mu_km = 398600.4418;  Re_km = 6378.1363;       % for printing only
DU = Re_km;  TU = sqrt(Re_km^3/mu_km);  VU = DU/TU;

outdir = fullfile(pwd, 'results', 'dc_equatorial_e065');
if ~exist(outdir, 'dir'), mkdir(outdir); end

C = [0.165 0.471 0.839     % blue
     0.922 0.408 0.204     % orange
     0.106 0.686 0.478     % green
     0.929 0.631 0.000];   % amber
ink = [0.10 0.10 0.10];  earthC = [0.85 0.85 0.85];

%% ---------------------------------------------------------------------
%  1. Kepler starting state (seed)
%  ---------------------------------------------------------------------
seed = dc_seed_perigee(alt_p_km, e_orb, Re_km, mu);
p    = seed.a*(1 - e_orb^2);                   % semi-latus rectum [DU]
fprintf('=== 1. ORBIT AND KEPLER SEED ====================================\n');
fprintf('Elements: perigee alt %.1f km, e = %.2f, i = 0, RAAN = 0, argp = 0, nu = 0\n', alt_p_km, e_orb);
fprintf('  r_p = %.12f DU = %.4f km\n', seed.rp, seed.rp_km);
fprintf('  r_a = %.12f DU = %.4f km  (apogee altitude %.1f km)\n', seed.ra, seed.ra_km, seed.ra_km - Re_km);
fprintf('  a   = %.12f DU = %.4f km\n', seed.a, seed.a_km);
fprintf('  Kepler period T = %.10f TU = %.4f min\n', seed.T, seed.T*TU/60);
fprintf('  Kepler perigee speed v_p = %.13f DU/TU = %.6f km/s\n\n', seed.vp, seed.vp*VU);

%% ---------------------------------------------------------------------
%  2. Differential correction
%  ---------------------------------------------------------------------
fprintf('=== 2. DIFFERENTIAL CORRECTION ==================================\n');
fprintf('Keep r_p fixed, correct v_p so the J2 orbit reaches r_a = %.4f km.\n', seed.ra_km);
dc = dc_correct_apogee(seed.rp, seed.ra, seed.vp, mu, J2, Re, opts);
fprintf(' iter   v_p [DU/TU]          change from Kepler [m/s]   apogee miss [km]\n');
for k = 1:size(dc.hist, 1)
    fprintf('  %d    %.13f     %+10.6f               %+.3e\n', dc.hist(k,1), dc.hist(k,2), ...
            (dc.hist(k,2) - seed.vp)*VU*1000, dc.hist(k,3)*DU);
end

% Independent check: at perigee and apogee rdot = 0, so energy and angular
% momentum h give one equation for h. Solving it needs no integration.
rp = seed.rp;  ra = seed.ra;
h2 = 2*((mu/rp - mu/ra) + J2*mu*Re^2/2*(1/rp^3 - 1/ra^3))/(1/rp^2 - 1/ra^2);
vp_alg = sqrt(h2)/rp;
fprintf('v_p by correction = %.15f\nv_p algebraic     = %.15f   diff %.1e\n', dc.vp, vp_alg, dc.vp - vp_alg);
fprintf('Extra perigee speed needed because of J2: %.4f m/s\n', (dc.vp - seed.vp)*VU*1000);
fprintf('Radial period T_r = %.10f TU = %.4f h   (Kepler T = %.4f h)\n', dc.Tr, dc.Tr*TU/3600, seed.T*TU/3600);
d_brouwer = 3*pi*J2*Re^2/p^2;                  % first-order, CITED (Brouwer 1959)
fprintf('Apsidal advance per radial period: %.8f deg  (first-order Brouwer %.8f deg, CITED)\n\n', ...
        rad2deg(dc.delta), rad2deg(d_brouwer));

%% ---------------------------------------------------------------------
%  3. One radial period: closure and monodromy
%  ---------------------------------------------------------------------
fprintf('=== 3. ONE RADIAL PERIOD =========================================\n');
X0 = dc.X0;
[~, X1, Phi] = dc_propagate(X0, [0 dc.Tr], mu, J2, Re, opts);
XT = X1(end, :).';
R6 = dc_rotz6(-dc.dtheta);
fprintf('|X(T_r) - X0|                  = %.3e DU (%.2f km)  -> not periodic in ECI\n', ...
        max(abs(XT - X0)), max(abs(XT(1:3) - X0(1:3)))*DU);
fprintf('|Rz(-dtheta) X(T_r) - X0|      = %.3e DU  -> periodic in the turning frame\n', max(abs(R6*XT - X0)));

Mrot   = R6*Phi;
lam_M  = eig(Mrot);
lam_P  = eig(Phi);
Jsym   = [zeros(3) eye(3); -eye(3) zeros(3)];
fprintf('symplecticity |M''JM - J| = %.2e,  |det M - 1| = %.2e\n', ...
        max(max(abs(Mrot.'*Jsym*Mrot - Jsym))), abs(det(Mrot) - 1));
fprintf('eig of Phi(T_r) alone (NOT a monodromy):  |lambda| = %s\n', sprintf('%.4f ', sort(abs(lam_P), 'descend')));
fprintf('eig of Rz(-dtheta)*Phi (monodromy):       |lambda| = %s\n', sprintf('%.12f ', sort(abs(lam_M), 'descend')));

ip = [1 2 4 5];  op = [3 6];                   
Mi = Mrot(ip, ip) - eye(4);
sv = svd(Mi);
fprintf('In-plane:  singular values of (M - I) = %s\n', sprintf('%.2e ', sv));
fprintf('           rank(M - I) = %d  (two 2x2 Jordan blocks at lambda = 1)\n', sum(sv > 1e-8*sv(1)));
fprintf('           ||(M - I)^2|| / ||M - I||^2 = %.2e\n', norm(Mi*Mi)/norm(Mi)^2);
lam_op   = eig(Mrot(op, op));
ang_op   = max(abs(angle(lam_op)));
fprintf('Out-of-plane pair: |lambda| = %.15f, angle = %.10f rad\n', abs(lam_op(1)), ang_op);
fprintf('   first-order prediction 6*pi*J2/p^2 = %.10f rad  (turning-frame rate minus node rate; uses Brouwer, CITED)\n', ...
        2*d_brouwer);
fprintf('   coupling in-plane <-> out-of-plane: %.1e (exactly decoupled at i = 0)\n\n', ...
        max(max(abs([Mrot(ip, op); Mrot(op, ip).']))));

%% ---------------------------------------------------------------------
%  4. Long run: N_REVS radial periods, sampled at every perigee
%  ---------------------------------------------------------------------
fprintf('=== 4. %d RADIAL PERIODS (%.1f days) ============================\n', N_REVS, N_REVS*dc.Tr*TU/86400);
tk = (0:N_REVS).' * dc.Tr;
[~, Xk] = dc_propagate(X0, tk, mu, J2, Re, opts);
rk   = sqrt(sum(Xk(:,1:3).^2, 2));
wk   = unwrap(atan2(Xk(:,2), Xk(:,1)));         % longitude of perigee (we are AT perigee)
revs = (0:N_REVS).';
Ek   = 0.5*sum(Xk(:,4:6).^2, 2) - mu./rk - J2*mu*Re^2./(2*rk.^3);   % z = 0 throughout
hk   = Xk(:,1).*Xk(:,5) - Xk(:,2).*Xk(:,4);
fprintf('Perigee radius stays at r_p to %.2e km over all %d revolutions\n', max(abs(rk - rp))*DU, N_REVS);
fprintf('Longitude of perigee after %d revs: %.5f deg  (correction predicts %.5f, Brouwer %.5f)\n', ...
        N_REVS, rad2deg(wk(end)), rad2deg(N_REVS*dc.delta), rad2deg(N_REVS*d_brouwer));
fprintf('Energy drift %.1e, angular momentum drift %.1e\n', max(abs(Ek - Ek(1))), max(abs(hk - hk(1))));
dXN = dc_rotz6(-N_REVS*dc.dtheta)*Xk(end,:).' - X0;
fprintf('Turning-frame closure after %d revs: position %.3e DU = %.3f m\n\n', N_REVS, ...
        norm(dXN(1:3)), norm(dXN(1:3))*DU*1000);

%% ---------------------------------------------------------------------
%  5. Figures
%  ---------------------------------------------------------------------
ang = linspace(0, 2*pi, 361);

% --- Figure 1: what the correction does ------------------------------
f1 = dc_light_figure([60 80 1400 480]);

% (a) seed orbit vs corrected orbit, one radial period each, under J2
axA = subplot(1, 3, 1); hold(axA, 'on');
[~, Xs] = dc_propagate(seed.X0, linspace(0, dc.Tr, 3000), mu, J2, Re, opts);
[~, Xc] = dc_propagate(X0,      linspace(0, dc.Tr, 3000), mu, J2, Re, opts);
fill(axA, cos(ang), sin(ang), earthC, 'EdgeColor', 'none', 'DisplayName', 'Earth');
plot(axA, Xs(:,1), Xs(:,2), '--', 'Color', C(2,:), 'LineWidth', 1.5, 'DisplayName', 'Kepler seed, under J2');
plot(axA, Xc(:,1), Xc(:,2), '-',  'Color', C(1,:), 'LineWidth', 1.5, 'DisplayName', 'corrected');
plot(axA, rp, 0, 'o', 'Color', ink, 'MarkerFaceColor', ink, 'HandleVisibility', 'off');
text(axA, rp + 0.2, 0.4, 'perigee (fixed)');
axis(axA, 'equal');
xlabel(axA, 'x  [DU]'); ylabel(axA, 'y  [DU]');
title(axA, 'One radial period (curves overlap at this scale)');
legend(axA, 'Location', 'southoutside');
dc_style_axes(axA);

% (b) apogee close-up for every Newton iterate
axB = subplot(1, 3, 2); hold(axB, 'on');
phi_a = pi + dc.delta/2;                       % direction of the target apogee
tw    = linspace(0.46, 0.54, 401)*dc.Tr;
nit   = size(dc.hist, 1);
for k = 1:nit
    [~, Xw] = dc_propagate([rp; 0; 0; 0; dc.hist(k,2); 0], [0 tw], mu, J2, Re, opts);
    Xw  = Xw(2:end, :);
    thw = mod(atan2(Xw(:,2), Xw(:,1)), 2*pi);
    rw  = sqrt(sum(Xw(:,1:3).^2, 2));
    cols = [C(2,:); C(4,:); C(3,:); C(1,:)];   % seed orange, then amber, green, blue
    col  = cols(min(k, 4), :);
    if k == 1, nm = 'iteration 0 (Kepler seed)'; else, nm = sprintf('iteration %d', k-1); end
    plot(axB, (thw - phi_a)*ra*DU, (rw - ra)*DU, 'Color', col, 'LineWidth', 1.8, 'DisplayName', nm);
end
plot(axB, [-3000 3000], [0 0], ':', 'Color', ink, 'LineWidth', 1.2, 'DisplayName', 'target apogee radius');
xlim(axB, [-2500 2500]);
xlabel(axB, 'along-track distance from target apogee  [km]');
ylabel(axB, 'r - r_a (target)  [km]');
title(axB, 'Close-up at apogee: iterations 1-3 overlap on the target');
legend(axB, 'Location', 'southoutside');
dc_style_axes(axB);

% (c) convergence
axC = subplot(1, 3, 3);
miss_km = abs(dc.hist(:,3))*DU;
semilogy(axC, dc.hist(:,1), max(miss_km, 1e-14), 'o-', 'Color', C(1,:), ...
         'LineWidth', 1.8, 'MarkerFaceColor', C(1,:), 'MarkerSize', 8);
for k = 1:nit
    text(axC, dc.hist(k,1) + 0.08, miss_km(k)*3, sprintf('%.1e km', miss_km(k)));
end
xlim(axC, [-0.3 nit - 0.3]); set(axC, 'XTick', 0:nit-1);
xlabel(axC, 'Newton iteration'); ylabel(axC, '|apogee miss|  [km]');
title(axC, sprintf('Convergence: extra %.3f m/s at perigee', (dc.vp - seed.vp)*VU*1000));
dc_style_axes(axC);
print(f1, fullfile(outdir, 'fig1_correction.png'), '-dpng', '-r200');

% --- Figure 2: the rosette over 400 revolutions ------------------------
fprintf('Figure 2: integrating 11 sample revolutions finely...\n');
kshow = 0:40:N_REVS;                           % revolutions to draw (0-based)
cmap  = [linspace(C(1,1), C(2,1), 256).', linspace(C(1,2), C(2,2), 256).', linspace(C(1,3), C(2,3), 256).'];
f2 = dc_light_figure([60 80 1300 600]);
axR = subplot(1, 2, 1); hold(axR, 'on');
axT = subplot(1, 2, 2); hold(axT, 'on');
fill(axR, cos(ang), sin(ang), earthC, 'EdgeColor', 'none');
fill(axT, cos(ang), sin(ang), earthC, 'EdgeColor', 'none');
for k = kshow
    [tf, Xf] = dc_propagate(Xk(k+1,:).', linspace(0, dc.Tr, 1500), mu, J2, Re, opts);
    col = cmap(1 + round(255*k/N_REVS), :);
    plot(axR, Xf(:,1), Xf(:,2), 'Color', col, 'LineWidth', 1.3);
    a  = -dc.delta/dc.Tr*(tk(k+1) + tf);       % undo the apsidal turning
    plot(axT, cos(a).*Xf(:,1) - sin(a).*Xf(:,2), sin(a).*Xf(:,1) + cos(a).*Xf(:,2), ...
         'Color', col, 'LineWidth', 1.3);
end
for ax = [axR axT]
    plot(ax, 0, 0, '+', 'Color', ink);
    axis(ax, 'equal'); xlim(ax, [-5.3 5.3]); ylim(ax, [-5.3 5.3]);   % apogee is 4.94 DU in any direction
    colormap(ax, cmap); clim_set = [1 N_REVS + 1];
    if exist('clim', 'file'), clim(ax, clim_set); else, caxis(ax, clim_set); end 
    cb = colorbar(ax); cb.Label.String = 'revolution'; cb.Color = ink;
end
xlabel(axR, 'x  [DU]'); ylabel(axR, 'y  [DU]');
title(axR, sprintf('ECI: every 40th revolution; perigee turns %.4f deg/rev', rad2deg(dc.delta)));
xlabel(axT, 'x''  [DU]'); ylabel(axT, 'y''  [DU]');
title(axT, 'Frame turning with the ellipse: all on one curve');
dc_style_axes(axR); dc_style_axes(axT);
print(f2, fullfile(outdir, 'fig2_rosette_400revs.png'), '-dpng', '-r200');

% --- Figure 3: apsidal drift --------------------------------------------
f3 = dc_light_figure([60 80 1300 480]);
axW = subplot(1, 2, 1); hold(axW, 'on');
plot(axW, revs, rad2deg(wk), '-', 'Color', C(1,:), 'LineWidth', 2.5, 'DisplayName', 'numerical (perigee position)');
plot(axW, revs, rad2deg(revs*d_brouwer), '--', 'Color', C(2,:), 'LineWidth', 1.5, 'DisplayName', 'first-order Brouwer (cited)');
xlabel(axW, 'revolution'); ylabel(axW, 'longitude of perigee  [deg]');
title(axW, sprintf('Perigee drift: %.3f deg after %d revs', rad2deg(wk(end)), N_REVS));
legend(axW, 'Location', 'northwest');
dc_style_axes(axW);
axE = subplot(1, 2, 2); hold(axE, 'on');
plot(axE, revs, rad2deg(wk - revs*dc.delta), '-', 'Color', C(1,:), 'LineWidth', 2, ...
     'DisplayName', 'numerical - correction prediction');
plot(axE, revs, rad2deg(revs*d_brouwer - revs*dc.delta), '--', 'Color', C(2,:), 'LineWidth', 1.5, ...
     'DisplayName', 'Brouwer - correction prediction');
xlabel(axE, 'revolution'); ylabel(axE, 'difference  [deg]');
title(axE, 'Differences (Brouwer misses the second-order J2^2 part)');
legend(axE, 'Location', 'southwest');
dc_style_axes(axE);
print(f3, fullfile(outdir, 'fig3_apsidal_drift.png'), '-dpng', '-r200');

% --- Figure 4: inside one radial period ---------------------------------
[t1, X1f] = dc_propagate(X0, linspace(0, dc.Tr, 4001), mu, J2, Re, opts);
r1  = sqrt(sum(X1f(:,1:3).^2, 2));
v1  = sqrt(sum(X1f(:,4:6).^2, 2));
aos = -mu./(2*(v1.^2/2 - mu./r1));             % osculating Kepler a
hv  = X1f(:,1).*X1f(:,5) - X1f(:,2).*X1f(:,4);
eos = sqrt(1 - hv.^2./(mu*aos));               % osculating Kepler e
th  = t1*TU/3600;
f4 = dc_light_figure([60 60 900 800]);
ax1 = subplot(3, 1, 1);
plot(ax1, th, (r1 - 1)*DU, 'Color', C(1,:), 'LineWidth', 2);
ylabel(ax1, 'altitude  [km]'); xlim(ax1, [0 th(end)]);
title(ax1, sprintf('One radial period (%.3f h): perigee %.1f km, apogee %.1f km', th(end), ...
      (rp - 1)*DU, (dc.ra - 1)*DU));
dc_style_axes(ax1);
ax2 = subplot(3, 1, 2);
plot(ax2, th, aos*DU, 'Color', C(2,:), 'LineWidth', 2);
ylabel(ax2, 'osculating a  [km]'); xlim(ax2, [0 th(end)]);
title(ax2, sprintf('Osculating semi-major axis varies by %.2f km within one orbit', (max(aos) - min(aos))*DU));
dc_style_axes(ax2);
ax3 = subplot(3, 1, 3);
plot(ax3, th, eos, 'Color', C(3,:), 'LineWidth', 2);
ylabel(ax3, 'osculating e'); xlabel(ax3, 'time since perigee  [h]'); xlim(ax3, [0 th(end)]);
title(ax3, sprintf('Osculating e: %.6f to %.6f (geometric e_{J2} = %.2f)', min(eos), max(eos), e_orb));
dc_style_axes(ax3);
print(f4, fullfile(outdir, 'fig4_one_period.png'), '-dpng', '-r200');

% --- Figure 5: monodromy eigenvalues -------------------------------------
f5 = dc_light_figure([60 80 1300 540]);
axM = subplot(1, 2, 1); hold(axM, 'on');
plot(axM, cos(ang), sin(ang), '-', 'Color', ink, 'LineWidth', 1, 'DisplayName', 'unit circle');
plot(axM, real(lam_P), imag(lam_P), 's', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', C(2,:), ...
     'MarkerFaceColor', C(2,:), 'DisplayName', '\Phi(T_r) alone (wrong frame)');
plot(axM, real(lam_M), imag(lam_M), 'o', 'MarkerSize', 12, 'LineWidth', 1.8, 'Color', C(1,:), ...
     'DisplayName', 'R_z(-\Delta\theta) \Phi(T_r)  (monodromy)');
axis(axM, 'equal'); xlim(axM, [-1.3 4.5]); ylim(axM, [-1.6 1.6]);
xlabel(axM, 'Re(\lambda)'); ylabel(axM, 'Im(\lambda)');
title(axM, sprintf('Without the rotation: fake |\\lambda| = %.2f', max(abs(lam_P))));
legend(axM, 'Location', 'southoutside');
dc_style_axes(axM);
axZ = subplot(1, 2, 2); hold(axZ, 'on');
tt = linspace(-0.01, 0.01, 400);
plot(axZ, cos(tt), sin(tt), '-', 'Color', ink, 'LineWidth', 1, 'DisplayName', 'unit circle');
plot(axZ, cos(2*d_brouwer)*[1 1], sin(2*d_brouwer)*[1 -1], 'x', 'MarkerSize', 14, 'LineWidth', 2, ...
     'Color', C(2,:), 'DisplayName', 'prediction 6\piJ_2/p^2');
plot(axZ, real(lam_op), imag(lam_op), 'o', 'MarkerSize', 10, 'LineWidth', 1.8, 'Color', C(1,:), ...
     'DisplayName', 'out-of-plane pair');
lam_in = eig(Mrot(ip, ip));
plot(axZ, real(lam_in), imag(lam_in), 'd', 'MarkerSize', 9, 'LineWidth', 1.5, 'Color', C(3,:), ...
     'MarkerFaceColor', C(3,:), 'DisplayName', 'in-plane (four at \lambda = 1)');
xlim(axZ, [0.99994 1.000008]); ylim(axZ, [-0.009 0.009]);
xlabel(axZ, 'Re(\lambda)'); ylabel(axZ, 'Im(\lambda)');
title(axZ, 'Monodromy, zoomed: all six on the unit circle');
legend(axZ, 'Location', 'west');
dc_style_axes(axZ);
print(f5, fullfile(outdir, 'fig5_monodromy.png'), '-dpng', '-r200');

%% ---------------------------------------------------------------------
%  6. Save
%  ---------------------------------------------------------------------
writetable(array2table([dc.hist(:,1), dc.hist(:,2), (dc.hist(:,2) - seed.vp)*VU*1000, dc.hist(:,3)*DU], ...
    'VariableNames', {'iteration', 'vp_DU_per_TU', 'dv_from_kepler_m_s', 'apogee_miss_km'}), ...
    fullfile(outdir, 'newton_history.csv'));
writetable(array2table([revs, tk, tk*TU/3600, rad2deg(wk), (rk - rp)*DU, Ek - Ek(1), hk - hk(1)], ...
    'VariableNames', {'rev', 't_TU', 't_hours', 'perigee_longitude_deg', 'perigee_radius_error_km', ...
    'energy_drift', 'h_drift'}), fullfile(outdir, 'per_revolution.csv'));
save(fullfile(outdir, 'dc_equatorial_e065.mat'), 'seed', 'dc', 'Phi', 'Mrot', 'lam_M', 'lam_P', ...
     'Xk', 'tk', 'wk', 'J2', 'N_REVS');
fprintf('Saved to %s\n', outdir);
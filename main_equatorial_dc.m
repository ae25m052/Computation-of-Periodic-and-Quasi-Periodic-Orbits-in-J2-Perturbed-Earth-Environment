%% main_equatorial_dc.m
%  Equatorial (i = 0) J2 orbits by differential correction, reproducing
%  Wang, Yuan, Zhao, Chen & Chen (2014), "Fourier Series Approximations to
%  J2-Bounded Equatorial Orbits", Math. Probl. Eng. 2014, 568318.
%
%  The paper uses closed-form elliptic integrals. Here everything comes
%  from integrating the equations of motion (j2_stm_eom.m, the validated
%  A matrix) and correcting the initial state with the STM.
%
%  CANONICAL UNITS THROUGHOUT:  DU = Re,  TU = sqrt(Re^3/mu),  mu = 1, Re = 1.
%  These are the same scales the paper uses (Sec. 3): r~ = r/Re,
%  t~ = t/sqrt(Re^3/mu), h~ = h/sqrt(Re*mu). All inputs below are already
%  canonical. km and s appear only in the printed output, for reading.
%
%  Sections
%    1. Units: paper vs this project
%    2. Critical circular orbit (a true periodic orbit) + Floquet multipliers
%    3. Pseudo-elliptic family, h = 1.1205, the five e values of Table 1
%    4. Largest e with r_min >= Re  (paper: e < 0.254)
%    5. Order-2 Fourier fit (Table 1) and order-1 deviation
%    6. Figures
%    7. Save results
%
%  Needs (same folder or on the path): j2_stm_eom.m, j2_A_matrix.m, and
%    eq_dc_circular.m   eq_dc_pseudo.m       eq_reference_apsides.m
%    eq_event_ycross.m  eq_event_apoapsis.m  eq_energy.m
%    eq_fourier_fit2.m  eq_fourier_eval.m    eq_sample_orbit.m
%    eq_light_figure.m  eq_style_axes.m

clear; clc; close all;

%% ---------------------------------------------------------------------
%  0. Settings
%  ---------------------------------------------------------------------
mu = 1;                 % canonical
Re = 1;                 % canonical
J2 = 1.08263e-3;        % THE PAPER'S VALUE (Sec. 3). The project value in
                        % j2_constants.m is 1.08262668e-3. The difference
                        % moves r0 by ~4e-9 DU (~0.03 m): 9th digit only.
h  = 1.1205;            % paper's fixed angular momentum (canonical)
rho_list = [0.80 0.75 0.70 0.65 0.60];     % r_min/r_max  ->  e = 0.1111 ... 0.25
paper_T1 = [ -9.99953e-2  9.99966e-2 -1.25461e-6  3.47383e-8
             -1.24991e-1  1.24994e-1 -2.99121e-6  1.11865e-7    % B0 printed as -1.24991e-2 (typo)
             -1.49988e-1  1.49991e-1 -3.18501e-6  2.36149e-7
             -1.74984e-1  1.74988e-1 -3.91574e-6  5.03017e-7
             -1.99789e-1  1.99984e-1 -4.71406e-6  9.65574e-7 ]; % paper Table 1: B0 B1 B2 minJ

tol  = 1e-13;
opts = odeset('RelTol', tol, 'AbsTol', tol);
opts_plot = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);   % long plotting runs only

% Dimensional scales -- used ONLY to print km / minutes alongside
mu_km = 398600.4418;  Re_km = 6378.1363;
DU = Re_km;  TU = sqrt(Re_km^3/mu_km);  VU = DU/TU;

outdir = fullfile(pwd, 'results', 'equatorial_dc');
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% ---------------------------------------------------------------------
%  1. Units
%  ---------------------------------------------------------------------
fprintf('=== 1. UNITS =====================================================\n');
fprintf('Paper (Sec. 3):  r~ = r/Re,  t~ = t/sqrt(Re^3/mu),  h~ = h/sqrt(Re*mu)\n');
fprintf('This project  :  DU = Re,    TU = sqrt(Re^3/mu),    h unit = DU^2/TU\n');
fprintf('  TU                = %.6f s\n', TU);
fprintf('  DU^2/TU           = %.4f km^2/s\n', DU^2/TU);
fprintf('  sqrt(Re*mu)       = %.4f km^2/s   (same number: identical scales)\n', sqrt(Re_km*mu_km));
fprintf('  energy unit VU^2  = %.6f km^2/s^2 = mu/Re = %.6f\n', VU^2, mu_km/Re_km);
fprintf('Note: the paper prints eps~ = eps/sqrt(Re/mu). That is dimensionally wrong;\n');
fprintf('its Figure 1 matches eps~ = eps/(mu/Re), which is our canonical energy.\n');
fprintf('h~ = %.4f  =  %.1f km^2/s\n\n', h, h*DU^2/TU);

%% ---------------------------------------------------------------------
%  2. Critical circular orbit
%  ---------------------------------------------------------------------
fprintf('=== 2. CRITICAL CIRCULAR ORBIT (h = %.4f) =========================\n', h);
[r0, Tc, histc] = eq_dc_circular(h, h^2/mu, mu, J2, Re, opts);   % Kepler guess r0 = h^2
fprintf(' iter        r0                 g = xdot at crossing     dg/dr0\n');
fprintf(' %3d   %.15f   %+.3e          %.6e\n', histc.');

r0_paper = (h^2 + sqrt(h^4 - 6*J2*mu^2*Re^2))/(2*mu);   % paper eq. (21)
T_paper  = 2*pi*r0_paper^2/h;                            % paper eq. (22)
fprintf('r0 (diff. correction) = %.15f DU = %.4f km (altitude %.2f km)\n', r0, r0*DU, (r0-1)*DU);
fprintf('r0 (paper eq. 21)     = %.15f DU    diff = %.2e\n', r0_paper, r0 - r0_paper);
fprintf('T  (diff. correction) = %.13f TU = %.4f min\n', Tc, Tc*TU/60);
fprintf('T  (paper eq. 22)     = %.13f TU    diff = %.2e\n', T_paper, Tc - T_paper);

% Monodromy matrix: a true periodic orbit, so Phi(T) IS the monodromy
% matrix and its eigenvalues ARE Floquet multipliers.
Xc0 = [r0; 0; 0; 0; h/r0; 0];
[~, Yc] = ode89(@(t,Y) j2_stm_eom(t, Y, mu, 1, J2, Re), [0 Tc], ...
                [Xc0; reshape(eye(6), 36, 1)], opts);
XcT = Yc(end, 1:6).';
M   = reshape(Yc(end, 7:42), 6, 6);
lam = eig(M);
Jsym = [zeros(3) eye(3); -eye(3) zeros(3)];
fprintf('closure |X(T) - X(0)|      = %.2e\n', max(abs(XcT - Xc0)));
fprintf('symplecticity |M''JM - J|  = %.2e,   |det M - 1| = %.2e\n', ...
        max(max(abs(M.'*Jsym*M - Jsym))), abs(det(M) - 1));
fprintf('Floquet multipliers:\n   Re(lambda)          Im(lambda)         |lambda|           angle [rad]\n');
fprintf('  %+.12f   %+.12f   %.15f   %+.10f\n', [real(lam) imag(lam) abs(lam) angle(lam)].');

chi  = sqrt(1 - 3*J2*mu*Re^2/(h^2*r0));                  % paper eq. (28), u0 = 1/r0
wzwt = sqrt((1 + 4.5*J2*Re^2/r0^2)/(1 + 1.5*J2*Re^2/r0^2));
fprintf('In-plane pair : predicted angle 2*pi*(1 - chi) = %.12f   (chi = %.12f, paper eq. 28)\n', ...
        2*pi*(1 - chi), chi);
fprintf('Out-of-plane  : predicted angle 2*pi*(wz/wth - 1) = %.12f  (DERIVED here, not in the paper)\n', ...
        2*pi*(wzwt - 1));
d2V = -2*mu/r0^3 - 6*J2*mu*Re^2/r0^5 + 3*h^2/r0^4;      % paper eq. (23)
fprintf('Paper eq. (23): d2Veff/dr2 at r0 = %.6f > 0  -> stable, consistent with |lambda| = 1\n\n', d2V);

%% ---------------------------------------------------------------------
%  3. Pseudo-elliptic family
%  ---------------------------------------------------------------------
fprintf('=== 3. PSEUDO-ELLIPTIC ORBITS (h = %.4f) ==========================\n', h);
nf  = numel(rho_list);
fam = struct([]);
for k = 1:nf
    rho = rho_list(k);
    e   = (1 - rho)/(1 + rho);
    d   = eq_dc_pseudo(h, e, mu, J2, Re, opts);
    ref = eq_reference_apsides(h, e, mu, J2, Re);

    % closure checks over one radial period
    [~, Y] = ode89(@(t,Y) j2_stm_eom(t, Y, mu, 1, J2, Re), [0 d.Tr], ...
                   [d.X0; reshape(eye(6), 36, 1)], opts);
    XT = Y(end, 1:6).';
    c = cos(d.dtheta);  s = sin(d.dtheta);
    Rz = [c -s 0; s c 0; 0 0 1];
    RX0 = [Rz*d.X0(1:3); Rz*d.X0(4:6)];

    fam(k).e = e;          fam(k).d = d;          fam(k).ref = ref;
    fam(k).rot_closure = max(abs(XT - RX0));
    fam(k).eci_closure = max(abs(XT - d.X0));

    fprintf('\ne_J2 = %.4f  (r_min/r_max = %.2f)   Newton |g| per iteration:', e, rho);
    fprintf(' %.1e', abs(d.hist(:,3)));
    fprintf('\n  r_min  = %.14f   ref %.14f   diff %+.1e   (%.3f km, alt %.1f km)\n', ...
            d.rmin, ref.rmin, d.rmin - ref.rmin, d.rmin*DU, (d.rmin-1)*DU);
    fprintf('  r_max  = %.14f   ref %.14f   diff %+.1e\n', d.rmax, ref.rmax, d.rmax - ref.rmax);
    fprintf('  energy = %.14f  ref %.14f   diff %+.1e\n', d.E, ref.E, d.E - ref.E);
    fprintf('  T_r    = %.13f   ref %.13f    diff %+.1e   (%.3f min)\n', ...
            d.Tr, ref.Tr, d.Tr - ref.Tr, d.Tr*TU/60);
    fprintf('  dtheta = %.13f   ref %.13f    diff %+.1e   (apsidal advance %.5f deg / radial period)\n', ...
            d.dtheta, ref.dtheta, d.dtheta - ref.dtheta, rad2deg(d.dtheta - 2*pi));
    fprintf('  |X(T_r) - Rz(dtheta) X0| = %.1e   (closes in the turning frame)\n', fam(k).rot_closure);
    fprintf('  |X(T_r) - X0|            = %.1e   (does NOT close in ECI)\n', fam(k).eci_closure);
end
fprintf('\n');

%% ---------------------------------------------------------------------
%  4. Largest e_J2 with the periapsis above the surface
%  ---------------------------------------------------------------------
fprintf('=== 4. e_J2 LIMIT (r_min = Re) ===================================\n');
X1 = [1; 0; 0; 0; h/1; 0];
[~, ~, te, ye] = ode89(@(t,Y) j2_stm_eom(t, Y, mu, 1, J2, Re), [0 50], ...
                       [X1; reshape(eye(6), 36, 1)], odeset(opts, 'Events', @eq_event_apoapsis));
rmax1 = norm(ye(end, 1:3));
e_lim = (rmax1 - 1)/(rmax1 + 1);
fprintf('Start at periapsis r = 1 DU, integrate to apoapsis: r_max = %.12f\n', rmax1);
fprintf('e_J2 limit = %.10f     paper: e_J2 < 0.254\n\n', e_lim);

%% ---------------------------------------------------------------------
%  5. Fourier series fits
%  ---------------------------------------------------------------------
fprintf('=== 5. FOURIER FITS (paper Sec. 5, Table 1) ======================\n');
nfit = 1001;      % samples over half a radial period, evenly spaced in time
fprintf('%d samples per half period. The paper does not state its number, and\n', nfit);
fprintf('min J grows with the number of samples, so compare B0, B1, B2 first.\n\n');
fprintf(' e_J2    |  B0 (ours / paper)          |  B1 (ours / paper)         |  B2 (ours / paper)          | min J (ours / paper)    | max dev ord.2 | max dev ord.1\n');
T1 = zeros(nf, 11);
for k = 1:nf
    d = fam(k).d;
    [~, th, r] = eq_sample_orbit(d.X0, d.Tr/2, nfit, mu, J2, Re, opts);
    [B, Jv, ra] = eq_fourier_fit2(th, r, d.rmin, d.rmax, d.dtheta/2);
    dev2 = max(abs(ra(:) - r)./r);
    B1o  = (d.rmax - d.rmin)/(2*d.rmax);                 % order 1, eq. (34)
    ra1  = eq_fourier_eval(th, [-B1o B1o], d.rmin, d.dtheta/2);
    dev1 = max(abs(ra1(:) - r)./r);
    P = paper_T1(k, :);
    fprintf(' %.4f  | %+.6e / %+.5e | %.6e / %.5e | %+.5e / %+.5e | %.2e / %.2e | %.1e       | %.1e\n', ...
            fam(k).e, B(1), P(1), B(2), P(2), B(3), P(3), Jv, P(4), dev2, dev1);
    T1(k, :) = [fam(k).e, B, Jv, P, dev2, dev1];
    fam(k).B = B;
end
fprintf('\ne_J2 implied by the paper''s own B1 (e = B1/(1 - B1)):');
fprintf(' %.6f', paper_T1(:,2)./(1 - paper_T1(:,2)));
fprintf('\nPaper Fig. 6 shows order-2 deviations up to ~2.2e-5; Fig. 5 shows order-1 up to ~0.067.\n\n');

%% ---------------------------------------------------------------------
%  6. Figures
%  ---------------------------------------------------------------------
% All figures: white background and explicit colours, so they look the
% same in MATLAB's light and dark themes (eq_light_figure, eq_style_axes).
C = [0.165 0.471 0.839     % blue
     0.922 0.408 0.204     % orange
     0.106 0.686 0.478     % green
     0.929 0.631 0.000];   % amber
ink    = [0.10 0.10 0.10];
earthC = [0.85 0.85 0.85];

% --- 6a. Paper Figure 1: energy vs initial apsis radius, fixed h -------
f1  = eq_light_figure([100 100 700 500]);
ax1 = axes(f1); hold(ax1, 'on');
r0v   = linspace(0.5, 1.5, 201);
hlist = [1.5172 1.3189 1.1205 0.9222];            % top curve first
for j = 1:numel(hlist)
    hh = hlist(j);
    Ev = arrayfun(@(rr) eq_energy([rr; 0; 0; 0; hh/rr; 0], mu, J2, Re), r0v);
    plot(ax1, r0v, Ev, 'Color', C(j,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('h = %.4f', hh));
end
plot(ax1, [0.5 1.5], [0 0], ':', 'Color', ink, 'LineWidth', 1.2, ...
     'DisplayName', '\epsilon = 0  (below: bounded, above: unbounded)');
xlabel(ax1, 'r_0  [DU]'); ylabel(ax1, '\epsilon  [canonical]');
xlim(ax1, [0.5 1.5]); ylim(ax1, [-1 3]);
title(ax1, 'Paper Fig. 1: energy vs initial apsis radius, fixed h');
legend(ax1, 'Location', 'northeast');
eq_style_axes(ax1);
print(f1, fullfile(outdir, 'fig1_energy_vs_r0.png'), '-dpng', '-r200');

% --- 6b. Paper Figure 8: r(t) for e = 0.25 with the order-2 fit --------
d = fam(end).d;  B = fam(end).B;
[t8, th8, r8] = eq_sample_orbit(d.X0, 50, 5001, mu, J2, Re, opts);
ra8  = eq_fourier_eval(th8, B, d.rmin, d.dtheta/2);
dev8 = abs(ra8(:) - r8)./r8;                       % the paper's Fig. 6 quantity

% Zoom window. The polar angle must be the TOTAL angle since t = 0
% (about 4.8 turns by t = 46.5), not the angle folded back into one turn:
% the Fourier series repeats every dtheta, which is slightly MORE than
% 2*pi, so folding would put r_app at the wrong point on the orbit. The
% whole turns are recovered from the continuous angle th8 above.
tz = [0, linspace(46.517, 46.519, 201)];
[~, Yz] = ode89(@(t,Y) j2_stm_eom(t, Y, mu, 1, J2, Re), tz, ...
                [d.X0; reshape(eye(6), 36, 1)], opts);
Yz = Yz(2:end, :);  tz = tz(2:end).';
rz    = sqrt(sum(Yz(:,1:3).^2, 2));
thw   = atan2(Yz(:,2), Yz(:,1));                   % folded into (-pi, pi]
turns = round((interp1(t8, th8, tz) - thw)/(2*pi));
thz   = thw + 2*pi*turns;                          % total angle
raz   = eq_fourier_eval(thz, B, d.rmin, d.dtheta/2);

f8  = eq_light_figure([100 100 1100 720]);
ax8 = subplot(2, 2, [1 2]); hold(ax8, 'on');
plot(ax8, t8, r8, 'Color', C(1,:), 'LineWidth', 2, 'DisplayName', 'numerical r(t)');
plot(ax8, t8, ra8, '--', 'Color', C(2,:), 'LineWidth', 1.5, 'DisplayName', 'order-2 Fourier r_{app}');
xlabel(ax8, 't  [TU]'); ylabel(ax8, 'r  [DU]'); xlim(ax8, [0 50]); ylim(ax8, [1 1.8]);
title(ax8, sprintf('Paper Fig. 8: e_{J2} = 0.25, h = %.4f', h));
legend(ax8, 'Location', 'northeast', 'Orientation', 'horizontal');
eq_style_axes(ax8);

axz = subplot(2, 2, 3); hold(axz, 'on');
plot(axz, tz, rz, 'Color', C(1,:), 'LineWidth', 2);
plot(axz, tz, raz, '--', 'Color', C(2,:), 'LineWidth', 1.5);
xlabel(axz, 't  [TU]'); ylabel(axz, 'r  [DU]'); xlim(axz, [46.517 46.519]);
set(axz, 'XTick', 46.517:0.0005:46.519);          % same ticks as the paper
xtickformat(axz, '%.4f');
title(axz, 'Zoom, t = 46.517 to 46.519 (curves overlap)');
eq_style_axes(axz);

axd = subplot(2, 2, 4);
semilogy(axd, t8, max(dev8, 1e-16), 'Color', C(3,:), 'LineWidth', 1.5);
xlabel(axd, 't  [TU]'); ylabel(axd, '|r_{app} - r| / r'); xlim(axd, [0 50]);
title(axd, sprintf('Fit error, max %.1e  (paper Fig. 6: ~2.2e-5)', max(dev8)));
eq_style_axes(axd);
print(f8, fullfile(outdir, 'fig8_radius_vs_time.png'), '-dpng', '-r200');
fprintf('Fig. 8 check: r(46.518) = %.6f DU (paper zoom shows ~1.3538)\n', interp1(tz, rz, 46.518));
fprintf('Fig. 8 zoom : max |r_app - r| in the window = %.1e DU\n', max(abs(raz(:) - rz)));

% --- 6c. The rosette in ECI, and the same orbit in the turning frame ----
delta = d.dtheta - 2*pi;                       % apsidal advance per radial period
kstep = round((pi/3)/delta);                   % revolutions for a 60 deg turn
nrev  = 2*kstep + 1;
nper  = 200;                                   % samples per radial period
fprintf('Rosette: %d radial periods (%.1f days). This run takes a minute or two...\n', ...
        nrev, nrev*d.Tr*TU/86400);
[tR, ~, ~, XR] = eq_sample_orbit(d.X0, nrev*d.Tr, nrev*nper + 1, mu, J2, Re, opts_plot);
ER = arrayfun(@(i) eq_energy(XR(i,:).', mu, J2, Re), (1:50:numel(tR)).');
fprintf('Energy drift over the rosette run: %.1e\n', max(abs(ER - ER(1))));

ang = linspace(0, 2*pi, 361);
fR  = eq_light_figure([100 100 1100 560]);
axR = subplot(1, 2, 1); hold(axR, 'on');
fill(axR, cos(ang), sin(ang), earthC, 'EdgeColor', 'none', 'DisplayName', 'Earth');
for j = 0:2
    kk  = j*kstep;                                          % revolution index (0-based)
    idx = (kk*nper + 1):((kk + 1)*nper + 1);
    plot(axR, XR(idx,1), XR(idx,2), 'Color', C(j+1,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('rev %d', kk + 1));
end
plot(axR, 0, 0, '+', 'Color', ink, 'HandleVisibility', 'off');
axis(axR, 'equal'); xlim(axR, [-1.9 1.9]); ylim(axR, [-1.9 1.9]);
xlabel(axR, 'x  [DU]'); ylabel(axR, 'y  [DU]');
title(axR, sprintf('ECI: the ellipse turns %.4f deg per revolution', rad2deg(delta)));
legend(axR, 'Location', 'southoutside', 'Orientation', 'horizontal');
eq_style_axes(axR);

axT = subplot(1, 2, 2); hold(axT, 'on');
a  = -delta/d.Tr * tR;                                     % rotate back at the apsidal rate
xr = cos(a).*XR(:,1) - sin(a).*XR(:,2);
yr = sin(a).*XR(:,1) + cos(a).*XR(:,2);
fill(axT, cos(ang), sin(ang), earthC, 'EdgeColor', 'none', 'DisplayName', 'Earth');
plot(axT, xr, yr, 'Color', C(1,:), 'LineWidth', 1.2, 'DisplayName', sprintf('all %d revs', nrev));
plot(axT, 0, 0, '+', 'Color', ink, 'HandleVisibility', 'off');
axis(axT, 'equal'); xlim(axT, [-1.9 1.9]); ylim(axT, [-1.9 1.9]);
xlabel(axT, 'x''  [DU]'); ylabel(axT, 'y''  [DU]');
title(axT, 'Turning frame: every revolution on one closed curve');
legend(axT, 'Location', 'southoutside', 'Orientation', 'horizontal');
eq_style_axes(axT);
print(fR, fullfile(outdir, 'rosette_eci_and_turning_frame.png'), '-dpng', '-r200');

% --- 6d. Floquet multipliers of the circular orbit ---------------------
% Sort the six multipliers into the three pairs by their angle.
ang_in  = 2*pi*(1 - chi);                   % paper eq. (28)
ang_oop = 2*pi*(wzwt - 1);                  % derived here
angl    = abs(angle(lam));
is_triv = angl < 1e-4;
is_in   = ~is_triv & abs(angl - ang_in) < abs(angl - ang_oop);
is_oop  = ~is_triv & ~is_in;

fF  = eq_light_figure([100 100 1100 500]);
axF = subplot(1, 2, 1); hold(axF, 'on');
tt  = linspace(-0.009, 0.009, 400);
plot(axF, cos(tt), sin(tt), '-', 'Color', ink, 'LineWidth', 1, 'DisplayName', 'unit circle');
plot(axF, real(lam(is_in)),  imag(lam(is_in)),  'o', 'MarkerSize', 10, 'LineWidth', 1.5, ...
     'Color', C(1,:), 'DisplayName', 'in-plane pair');
plot(axF, real(lam(is_oop)), imag(lam(is_oop)), 's', 'MarkerSize', 6, 'LineWidth', 1.5, ...
     'Color', C(2,:), 'MarkerFaceColor', C(2,:), 'DisplayName', 'out-of-plane pair');
plot(axF, real(lam(is_triv)), imag(lam(is_triv)), 'd', 'MarkerSize', 8, 'LineWidth', 1.5, ...
     'Color', C(3,:), 'MarkerFaceColor', C(3,:), 'DisplayName', 'trivial pair (\lambda = 1)');
xlim(axF, [0.99996 1.000006]); ylim(axF, [-0.009 0.009]);
xlabel(axF, 'Re(\lambda)'); ylabel(axF, 'Im(\lambda)');
title(axF, 'All six multipliers lie on the unit circle');
legend(axF, 'Location', 'west');
eq_style_axes(axF);

% Right panel: the upper in-plane and out-of-plane multipliers overlap on
% the left, so plot each one's offset from its predicted angle instead.
axG = subplot(1, 2, 2); hold(axG, 'on');
up_in  = lam(is_in  & imag(lam) > 0);
up_oop = lam(is_oop & imag(lam) > 0);
plot(axG, [0 0], [-1e3 1e3], '--', 'Color', C(1,:), 'LineWidth', 1, 'HandleVisibility', 'off');
plot(axG, (ang_oop - ang_in)*1e6*[1 1], [-1e3 1e3], '--', 'Color', C(2,:), 'LineWidth', 1, 'HandleVisibility', 'off');
plot(axG, (angle(up_in)  - ang_in)*1e6, (abs(up_in)  - 1)*1e12, 'o', 'MarkerSize', 10, ...
     'LineWidth', 1.5, 'Color', C(1,:), 'DisplayName', 'in-plane (predicted: blue dashed)');
plot(axG, (angle(up_oop) - ang_in)*1e6, (abs(up_oop) - 1)*1e12, 's', 'MarkerSize', 8, ...
     'LineWidth', 1.5, 'Color', C(2,:), 'MarkerFaceColor', C(2,:), ...
     'DisplayName', 'out-of-plane (predicted: orange dashed)');
yl = max(1, 1.3*max(abs([abs(up_in); abs(up_oop)] - 1))*1e12);
xlim(axG, [-8 2]); ylim(axG, [-yl yl]);
xlabel(axG, 'angle - 2\pi(1-\chi)   [10^{-6} rad]');
ylabel(axG, '|\lambda| - 1   [10^{-12}]');
title(axG, 'Upper pair, separated: measured vs predicted');
legend(axG, 'Location', 'southoutside');
eq_style_axes(axG);
print(fF, fullfile(outdir, 'floquet_circular.png'), '-dpng', '-r200');

%% ---------------------------------------------------------------------
%  7. Save
%  ---------------------------------------------------------------------
famtab = zeros(nf, 12);
for k = 1:nf
    famtab(k,:) = [fam(k).e, fam(k).d.rmin, fam(k).d.rmax, fam(k).d.E, fam(k).d.Tr, ...
                   fam(k).d.dtheta, fam(k).ref.rmin, fam(k).ref.Tr, fam(k).ref.dtheta, ...
                   fam(k).rot_closure, fam(k).eci_closure, size(fam(k).d.hist, 1) - 1];
end
writetable(array2table(famtab, 'VariableNames', {'e_J2','rmin','rmax','energy','Tr', ...
    'dtheta','rmin_ref','Tr_ref','dtheta_ref','rot_closure','eci_closure','newton_steps'}), ...
    fullfile(outdir, 'pseudo_elliptic_family.csv'));
writetable(array2table(T1, 'VariableNames', {'e_J2','B0','B1','B2','minJ', ...
    'B0_paper','B1_paper','B2_paper','minJ_paper','maxdev_order2','maxdev_order1'}), ...
    fullfile(outdir, 'fourier_table1.csv'));
save(fullfile(outdir, 'equatorial_dc.mat'), 'r0', 'Tc', 'M', 'lam', 'fam', 'T1', 'e_lim', 'J2', 'h');
fprintf('\nSaved to %s\n', outdir);
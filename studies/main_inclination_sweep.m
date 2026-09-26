%% MAIN_INCLINATION_SWEEP
%  Sweep the inclination from 0 to 180 degrees. At each inclination:
%
%     (a) propagate for N_REVS revolutions at epsilon = 0 (pure Kepler) and
%         at epsilon = 1 (full J2), and plot each in its own folder;
%     (b) test closure, and measure the secular drift of the node and the
%         apsides, which is where the closure question is actually settled;
%     (c) propagate the state transition matrix over ONE KEPLERIAN PERIOD
%         for both models and take its eigenvalues.
%
%  Every routine this script needs lives in its own function file. There
%  are no local functions here -- the script only sets up, loops, plots and
%  saves.
%
%  ---------------------------------------------------------------------
%  REQUIRES on the MATLAB path
%
%  New function files (supplied alongside this script):
%     j2_constants.m            seed_orbit_from_perigee.m
%     seed_initial_state.m      prop_state_sol.m
%     prop_stm.m                rv2coe_safe.m
%     element_history.m         stm_report.m
%     brouwer_rates.m           inc_tag.m
%     thin_index.m              save_figure_png.m
%     format_orbit_axes.m
%
%  Existing project files (used unchanged):
%     coe2rv.m        twobody_eom.m          j2_A_matrix.m
%     j2_stm_eom.m    plot_earth.m           plot_orbit_gradient.m
%     style_light.m   export_for_slides.m
%
%  DELIBERATELY NOT USED: rv2coe.m, which returns NaN for argp at i = 0
%  because it divides by the norm of the zero node vector. rv2coe_safe.m is
%  the guarded replacement. rv2coe.m itself is untouched and still correct
%  everywhere away from the two equatorial endpoints.
%  ---------------------------------------------------------------------
%
%  M.Tech. project, Raj Khismatrao (AE25M052).
%  Guide: Prof. Joel George Manathara, IIT Madras.

clear; clc; close all;

%% =====================================================================
%  1. CONFIGURATION   <<<< everything you might want to change is here
%  =====================================================================

% ---- Run the short version first ------------------------------------
%  Leave this TRUE the first time. It uses 3 inclinations and 20
%  revolutions and finishes quickly, which is enough to confirm that every
%  folder, file name and figure comes out right. Only then set it FALSE for
%  the full run.
%
%  Setting it FALSE re-runs EVERYTHING from scratch, including the three
%  inclinations the quick test already covered -- this time at the full
%  revolution count. The quick-test PNGs are overwritten, which is what you
%  want: 20-revolution figures are throwaway.
QUICK_TEST = false;

% ---- Where the output goes -------------------------------------------
OUT_ROOT = 'D:\MATLAB\ECI_project\results\inclination_sweep';

% ---- The orbit (as specified by the guide) ---------------------------
HP_KM    = 300;      % perigee ALTITUDE above the surface [km]
ECC      = 0.65;     % eccentricity -- "more elliptic"
RAAN_DEG = 0;        % the guide set both of these to zero
ARGP_DEG = 0;
NU_DEG   = 0;        % "any value"; perigee is the natural choice

% ---- How long, and at what accuracy ----------------------------------
N_REVS   = 400;      % revolutions for the orbital-variation study
TOL_TRAJ = 1e-11;    % tolerance for the long trajectory runs
TOL_STM  = 1e-13;    % tolerance for the one-period STM runs

% ---- Which inclinations ----------------------------------------------
%  Two separate lists, on purpose.
%
%  INC_LIST drives the expensive part: long propagations and three figures
%  each. One-degree steps -- 181 values, plus the two critical inclinations
%  inserted, so 183 in all.
%
%  INC_EIG drives the cheap part: one Keplerian period each, no figures,
%  eigenvalues only. Two-and-a-half-degree steps -- 73 values plus the two
%  critical ones, so 75 in all. Finer, so the eigenvalue-versus-inclination
%  curve comes out smooth instead of joining fifteen-degree dots.
%
%  Both lists MUST contain the critical inclination. On a one-degree grid
%  the nearest points are 63 and 64 deg, and 63.434949 falls between them,
%  so it would otherwise be missed -- and it is the one inclination the
%  guide asked about specifically.
i_crit   = acosd(1/sqrt(5));                        % 63.434949 deg
INC_LIST = sort([0:1:180, i_crit, 180 - i_crit]);
INC_EIG  = sort([0:2.5:180, i_crit, 180 - i_crit]);

% ---- Plotting --------------------------------------------------------
MAX_PLOT_PTS = 12000;   % points kept for a 3D trace (see thin_index.m)
SHOW_FIGURES = false;   % false = build figures off-screen. Keep it false:
                        % a full run makes 549 of them.

if QUICK_TEST
    N_REVS   = 20;
    INC_LIST = [0, 45, i_crit];
    INC_EIG  = sort([0:15:180, i_crit, 180 - i_crit]);
    fprintf('*** QUICK_TEST is ON: %d revolutions, %d inclinations ***\n\n', ...
            N_REVS, numel(INC_LIST));
end

%% =====================================================================
%  2. SETUP
%  =====================================================================

needed = {'j2_constants','seed_orbit_from_perigee','seed_initial_state', ...
          'prop_state_sol','prop_stm','rv2coe_safe','element_history', ...
          'stm_report','brouwer_rates','inc_tag','thin_index', ...
          'save_figure_png','format_orbit_axes', ...
          'coe2rv','twobody_eom','j2_A_matrix','j2_stm_eom', ...
          'plot_earth','plot_orbit_gradient','style_light','export_for_slides'};
missing = needed(cellfun(@(f) isempty(which(f)), needed));
if ~isempty(missing)
    error(['These files are not on the MATLAB path:\n   %s\n\n' ...
           'If you have just copied them in, MATLAB cannot see new files ' ...
           'until it is RESTARTED.'], strjoin(missing, ', '));
end

c = j2_constants();
s = seed_orbit_from_perigee(HP_KM, ECC);

vis = 'off';
if SHOW_FIGURES
    vis = 'on';
end

dirK = fullfile(OUT_ROOT, 'kepler');
dirJ = fullfile(OUT_ROOT, 'j2');
dirV = fullfile(OUT_ROOT, 'eigenvalues');
dirD = fullfile(OUT_ROOT, 'data');

% MATLAB .fig copies of the J2 trajectory figures only, so those can be
% reopened and rotated. The Kepler and eigenvalue figures stay PNG-only.
dirFJ = fullfile(OUT_ROOT, 'j2_fig');

all_dirs = {OUT_ROOT, dirK, dirJ, dirV, dirD, dirFJ};
for kd = 1:numel(all_dirs)
    if ~exist(all_dirs{kd}, 'dir')
        mkdir(all_dirs{kd});
    end
end

fprintf('=========================================================\n');
fprintf(' INCLINATION SWEEP\n');
fprintf('=========================================================\n');
fprintf(' perigee altitude  : %.1f km   (r_p = %.4f km)\n', s.hp_km, s.rp_km);
fprintf(' eccentricity      : %.4f\n', s.ecc);
fprintf(' semi-major axis   : %.4f km  (follows from h_p and e)\n', s.a_km);
fprintf(' apogee altitude   : %.1f km\n', s.ha_km);
fprintf(' Keplerian period  : %.4f min  (%.6f TU)\n', s.T_min, s.T);
fprintf(' %d revolutions    : %.2f days\n', N_REVS, N_REVS*s.T*c.TU/86400);
fprintf(' critical incl.    : %.6f deg  and  %.6f deg\n', ...
        c.i_crit_deg, c.i_crit_retro);
fprintf(' inclinations      : %d with figures, %d for the eigenvalue curve\n', ...
        numel(INC_LIST), numel(INC_EIG));
fprintf(' output root       : %s\n', OUT_ROOT);
fprintf('=========================================================\n\n');

%% =====================================================================
%  3. PHASE A -- fine eigenvalue sweep (one Keplerian period, no figures)
%  =====================================================================

fprintf('PHASE A: eigenvalues over one Keplerian period, %d inclinations\n', ...
        numel(INC_EIG));

nE       = numel(INC_EIG);
maxabs_K = zeros(nE,1);   maxabs_J = zeros(nE,1);
symp_K   = zeros(nE,1);   symp_J   = zeros(nE,1);
absev_K  = zeros(nE,6);   absev_J  = zeros(nE,6);
pair_J   = zeros(nE,1);
jres_K   = zeros(nE,1);   rank_K   = zeros(nE,1);

tA = tic;
for kk = 1:nE

    X0 = seed_initial_state(s, INC_EIG(kk), RAAN_DEG, ARGP_DEG, NU_DEG);

    [~, PhiK] = prop_stm(X0, s.T, 0, TOL_STM);
    [~, PhiJ] = prop_stm(X0, s.T, 1, TOL_STM);

    RK = stm_report(PhiK);
    RJ = stm_report(PhiJ);

    maxabs_K(kk)  = RK.maxabs;    maxabs_J(kk)  = RJ.maxabs;
    symp_K(kk)    = RK.symp;      symp_J(kk)    = RJ.symp;
    absev_K(kk,:) = RK.absev.';   absev_J(kk,:) = RJ.absev.';
    pair_J(kk)    = RJ.pair_err;
    jres_K(kk)    = RK.jordan_resid;
    rank_K(kk)    = RK.rank_PhiI;

    if mod(kk, 10) == 0 || kk == nE
        fprintf('   %3d/%3d   i = %7.3f deg   max|lambda|_J2 = %10.6f\n', ...
                kk, nE, INC_EIG(kk), RJ.maxabs);
    end
end
fprintf('PHASE A done in %.1f s\n\n', toc(tA));

% ---- CHECK 1: mirror symmetry ---------------------------------------
%  With RAAN = argp = 0 the dynamics depend on inclination only through
%  cos^2(i), so the eigenvalue curve must be symmetric about 90 deg. The
%  achievable level is set by the integrator error, i.e. by the
%  symplecticity figure below -- not by machine precision.
fprintf('CHECK 1  mirror symmetry of max|lambda| about 90 deg : %.3e\n', ...
        max(abs(maxabs_J - flipud(maxabs_J))));
fprintf('         (expect the same order as the symplecticity, ~1e-10)\n');
fprintf('CHECK 1  worst symplecticity over the sweep          : %.3e\n\n', ...
        max(symp_J));

% ---- CHECK 2: the Keplerian STM is a shear --------------------------
%  Tested STRUCTURALLY, not through the eigenvalues. See stm_report.m:
%  computing the eigenvalues of a defective matrix is ill-conditioned by
%  nature, so a Keplerian STM that is perfectly correct can still report
%  eigenvalues that miss 1 by ~1e-6. rank(Phi-I) = 1 and (Phi-I)^2 = 0 are
%  the trustworthy tests, and they go nowhere near that difficulty.
fprintf('CHECK 2  rank(Phi - I) at eps = 0, over the sweep    : min %d, max %d\n', ...
        min(rank_K), max(rank_K));
fprintf('         (must be 1 everywhere: a single Jordan block)\n');
fprintf('CHECK 2  ||(Phi-I)^2|| / ||Phi-I||^2 at eps = 0      : %.3e\n', ...
        max(jres_K));
fprintf('         (must be small: the shear squares to zero)\n');
fprintf('CHECK 2  max | |lambda| - 1 | at eps = 0             : %.3e\n', ...
        max(abs(maxabs_K - 1)));
fprintf('         (NOT a failure if this is ~1e-6; it is the defective\n');
fprintf('          eigenvalue responding to rounding like sqrt of it)\n\n');

% ---- CHECK 3: reciprocal pairing at eps = 1 -------------------------
fprintf('CHECK 3  worst reciprocal-pair error at eps = 1      : %.3e\n', ...
        max(pair_J));
fprintf('         (the structural check that survives into the perturbed\n');
fprintf('          problem, where the eigenvalues are no longer all 1)\n\n');

% index of the polar case, used as the reference for the normalised plots
[~, i90] = min(abs(INC_EIG - 90));

%% =====================================================================
%  4. PHASE B -- per-inclination trajectories, closure and figures
%  =====================================================================

fprintf('PHASE B: %d revolutions at %d inclinations, 3 figures each\n', ...
        N_REVS, numel(INC_LIST));

nI    = numel(INC_LIST);
t_rev = (0:N_REVS) * s.T;          % sample once per Keplerian period
t_end = N_REVS * s.T;

res = struct('i_deg',{},'d_argp',{},'d_RAAN',{},'d_argp_pred',{}, ...
             'd_RAAN_pred',{},'clos_end',{},'clos_min',{},'clos_min_rev',{}, ...
             'maxabs_K',{},'maxabs_J',{},'symp_J',{},'equatorial',{});

tB = tic;
for kk = 1:nI

    inc = INC_LIST(kk);
    tag = inc_tag(inc);
    fprintf('\n--- i = %8.4f deg   (tag %s)   [%d/%d] ---\n', inc, tag, kk, nI);

    X0 = seed_initial_state(s, inc, RAAN_DEG, ARGP_DEG, NU_DEG);

    % ---- propagate both models once, then sample them several ways ----
    solK = prop_state_sol(X0, t_end, 0, TOL_TRAJ);
    solJ = prop_state_sol(X0, t_end, 1, TOL_TRAJ);

    XK = solK.y;   tK = solK.x;     % the integrator's own steps
    XJ = solJ.y;   tJ = solJ.x;

    XJrev = deval(solJ, t_rev);     % once per period, for the elements

    fprintf('    integrator steps: Kepler %d, J2 %d\n', numel(tK), numel(tJ));

    % ---- secular drift -------------------------------------------------
    HJ = element_history(XJrev, c.mu);
    SR = brouwer_rates(s, inc);

    d_argp_pred = rad2deg(SR.dargp) * t_end;
    d_RAAN_pred = rad2deg(SR.dRAAN) * t_end;

    if HJ.equatorial
        % An equatorial orbit has no ascending node. rv2coe_safe therefore
        % reports the LONGITUDE OF PERIGEE and sets RAAN to zero, so the
        % first-order prediction has to be the matching COMBINATION of the
        % two Brouwer rates rather than the argp rate on its own:
        %     prograde   varpi = RAAN + argp
        %     retrograde varpi = RAAN - argp   (motion runs the other way)
        % Comparing the measured longitude of perigee against the argp rate
        % alone is comparing two different angles, and disagrees by roughly
        % a factor of two.
        angname = 'longitude of perigee';
        if inc < 90
            d_argp_pred = d_argp_pred + d_RAAN_pred;
        else
            d_argp_pred = d_argp_pred - d_RAAN_pred;
        end
        d_RAAN_pred = 0;
        fprintf('    EQUATORIAL: no ascending node exists, so the angle below\n');
        fprintf('                is the longitude of perigee, and the Brouwer\n');
        fprintf('                figure is the matching combination of rates.\n');
    else
        angname = 'argument of perigee';
    end

    fprintf('    d(%-21s) : %+10.4f deg   (Brouwer 1st order %+10.4f)\n', ...
            angname, HJ.d_argp_deg, d_argp_pred);
    fprintf('    d(%-21s) : %+10.4f deg   (Brouwer 1st order %+10.4f)\n', ...
            'RAAN', HJ.d_RAAN_deg, d_RAAN_pred);

    % ---- closure -------------------------------------------------------
    %  How close does the J2 orbit ever come to its own starting point,
    %  sampled at whole Keplerian periods? Closure in its bluntest form.
    dclose = vecnorm(XJrev(1:3,2:end) - X0(1:3)) * c.DU;
    [clos_min, imin] = min(dclose);
    clos_end = dclose(end);
    fprintf('    closure: final %.1f km, best %.1f km at revolution %d\n', ...
            clos_end, clos_min, imin);

    % ---- STM over one Keplerian period ---------------------------------
    [~, PhiK] = prop_stm(X0, s.T, 0, TOL_STM);
    [~, PhiJ] = prop_stm(X0, s.T, 1, TOL_STM);
    RK = stm_report(PhiK);
    RJ = stm_report(PhiJ);
    fprintf('    max|lambda|: Kepler %.10f, J2 %.6f   (symplecticity %.2e)\n', ...
            RK.maxabs, RJ.maxabs, RJ.symp);

    % =================================================================
    %  FIGURE 1 -- Keplerian trajectory
    % =================================================================
    iK = thin_index(numel(tK), MAX_PLOT_PTS);

    f1 = figure('Visible', vis, 'Color', 'w', 'Name', ['kepler ' tag]);
    plot_earth(c.Re); hold on;
    plot3(XK(1,iK), XK(2,iK), XK(3,iK), '-', ...
          'Color', [0.00 0.35 0.85], 'LineWidth', 1.6);
    format_orbit_axes(XK, sprintf('Kepler (\\epsilon = 0),  i = %.2f^\\circ,  %d revs', ...
                      inc, N_REVS));
    style_light(f1, 14);
    save_figure_png(f1, dirK, ['kepler_inc' tag], struct('Width',9,'Height',8));
    close(f1);

    % =================================================================
    %  FIGURE 2 -- J2 trajectory, coloured early -> late
    % =================================================================
    iJ = thin_index(numel(tJ), MAX_PLOT_PTS);

    f2 = figure('Visible', vis, 'Color', 'w', 'Name', ['j2 ' tag]);
    plot_earth(c.Re); hold on;
    plot_orbit_gradient(XJ(:,iJ), tJ(iJ)*c.TU/86400, 'turbo', 1.2);
    % The first revolution is picked out in black, so the precession of the
    % ellipse away from where it began is visible at a glance. At the
    % critical inclination the swept surface is a clean fan about the polar
    % axis; elsewhere the ellipse also tumbles inside its own plane and the
    % trace smears into a thick shell.
    i1 = tJ <= s.T;
    hfirst = plot3(XJ(1,i1), XJ(2,i1), XJ(3,i1), 'k-', 'LineWidth', 2.2);
    cb = colorbar;
    cb.Label.String = 'time  [days]';
    format_orbit_axes(XJ, sprintf('J_2 (\\epsilon = 1),  i = %.2f^\\circ,  %d revs', ...
                      inc, N_REVS));
    legend(hfirst, {'first revolution'}, 'Location', 'northeast');
    style_light(f2, 14);
    save_figure_png(f2, dirJ, ['j2_inc' tag], struct('Width',10,'Height',8), dirFJ);
    close(f2);

    % =================================================================
    %  FIGURE 3 -- eigenvalues of Phi(T), three views
    %
    %  The third panel exists because the first two look almost identical
    %  from one inclination to the next: max|lambda| only moves from 90.73
    %  at equatorial to 87.51 at polar, about 3.5 per cent, which on a log
    %  axis spanning 0.011 to 90 is roughly one pixel. Panel 3 places this
    %  inclination on the sweep curve with a linear, zoomed axis, so the
    %  variation is actually visible and each figure is distinguishable.
    % =================================================================
    f3 = figure('Visible', vis, 'Color', 'w', 'Name', ['eig ' tag]);

    subplot(1,3,1);
    th = linspace(0, 2*pi, 400);
    plot(cos(th), sin(th), 'k-', 'LineWidth', 1.0); hold on;
    hk = plot(real(RK.ev), imag(RK.ev), 'o', 'MarkerSize', 11, ...
              'Color', [0.00 0.35 0.85], 'LineWidth', 1.8);
    hj = plot(real(RJ.ev), imag(RJ.ev), 'x', 'MarkerSize', 12, ...
              'Color', [0.85 0.10 0.10], 'LineWidth', 2.0);
    grid on; axis equal;
    xlim([-1.6 1.6]); ylim([-1.3 1.3]);
    xlabel('Re \lambda'); ylabel('Im \lambda');
    legend([hk hj], {'Kepler', 'J_2'}, 'Location', 'south');
    title('near the unit circle');

    subplot(1,3,2);
    semilogy(1:6, RK.absev, 'o-', 'Color', [0.00 0.35 0.85], ...
             'LineWidth', 1.8, 'MarkerSize', 9); hold on;
    semilogy(1:6, RJ.absev, 'x-', 'Color', [0.85 0.10 0.10], ...
             'LineWidth', 2.0, 'MarkerSize', 11);
    yline(1, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    % Print the two extreme magnitudes on the figure, so it carries the
    % numbers and can be read without opening the CSV.
    text(6, RJ.absev(6), sprintf('  %.3f', RJ.absev(6)), ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    text(1, RJ.absev(1), sprintf('%.5f  ', RJ.absev(1)), ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    grid on; xlim([0.5 6.5]); xticks(1:6); ylim([5e-3 3e2]);
    xlabel('eigenvalue index'); ylabel('|\lambda|');
    legend('Kepler', 'J_2', 'Location', 'northwest');
    title('magnitudes, log scale');

    subplot(1,3,3);
    plot(INC_EIG, maxabs_J, '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 2.0);
    hold on;
    plot(inc, RJ.maxabs, 'ko', 'MarkerSize', 13, 'LineWidth', 2.5, ...
         'MarkerFaceColor', [1 1 0.3]);
    xline(c.i_crit_deg, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xline(c.i_crit_retro, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    grid on; xlim([-10 190]); xticks(0:45:180);
    xlabel('inclination  [deg]'); ylabel('max |\lambda|');
    title('position on the sweep');

    sgt = sgtitle(sprintf(['Eigenvalues of \\Phi(T) over ONE Keplerian period,  ' ...
        'i = %.2f^\\circ'], inc));
    set(sgt, 'FontWeight', 'bold');
    style_light(f3, 12);
    save_figure_png(f3, dirV, ['eig_inc' tag], ...
                    struct('Width',17,'Height',6.5,'FontSize',14,'TitleSize',15));
    close(f3);

    % ---- record --------------------------------------------------------
    res(kk).i_deg        = inc;
    res(kk).d_argp       = HJ.d_argp_deg;
    res(kk).d_RAAN       = HJ.d_RAAN_deg;
    res(kk).d_argp_pred  = d_argp_pred;
    res(kk).d_RAAN_pred  = d_RAAN_pred;
    res(kk).clos_end     = clos_end;
    res(kk).clos_min     = clos_min;
    res(kk).clos_min_rev = imin;
    res(kk).maxabs_K     = RK.maxabs;
    res(kk).maxabs_J     = RJ.maxabs;
    res(kk).symp_J       = RJ.symp;
    res(kk).equatorial   = HJ.equatorial;

end
fprintf('\nPHASE B done in %.1f s\n\n', toc(tB));

%% =====================================================================
%  5. SUMMARY FIGURES
%  =====================================================================

iv = [res.i_deg];

% ---- max |lambda| versus inclination ---------------------------------
%  Two panels. The left one is the honest overview and shows that the
%  Keplerian curve sits at 1 while the J2 curve sits near 88. The right one
%  is the same J2 curve alone on a LINEAR, zoomed axis -- which is the only
%  way the real structure (a shallow U with its minimum at polar) becomes
%  visible at all. On the left panel that whole structure is one pixel
%  thick.
fs1 = figure('Visible', vis, 'Color', 'w', 'Name', 'eig sweep maxabs');

subplot(1,2,1);
hJ = semilogy(INC_EIG, maxabs_J, '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 2.2);
hold on;
hK = semilogy(INC_EIG, maxabs_K, '-', 'Color', [0.00 0.35 0.85], 'LineWidth', 2.2);
hC = xline(c.i_crit_deg, 'k--', 'LineWidth', 1.4);
% The remaining reference lines are hidden from the legend, otherwise each
% one claims an entry and pushes the real labels onto the wrong curves.
h2 = xline(c.i_crit_retro, 'k--', 'LineWidth', 1.4);
set(h2, 'HandleVisibility', 'off');
grid on; xlim([-10 190]); xticks(0:45:180); ylim([0.5 200]);
xlabel('inclination  [deg]'); ylabel('max |\lambda|  of  \Phi(T)');
legend([hJ hK hC], {'J_2 (\epsilon = 1)', 'Kepler (\epsilon = 0)', ...
       'critical incl.'}, 'Location', 'east');
title('overview, log scale');

subplot(1,2,2);
plot(INC_EIG, maxabs_J, '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 2.4);
hold on;
xline(c.i_crit_deg, 'k--', 'LineWidth', 1.4, 'HandleVisibility', 'off');
xline(c.i_crit_retro, 'k--', 'LineWidth', 1.4, 'HandleVisibility', 'off');
xline(90, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
grid on; xlim([-10 190]); xticks(0:45:180);
xlabel('inclination  [deg]'); ylabel('max |\lambda|  of  \Phi(T)');
title('J_2 alone, linear and zoomed');

sgt = sgtitle({'Largest eigenvalue of \Phi(T) over one Keplerian period', ...
       'smooth, symmetric about 90^\circ, and with NO feature at the critical inclination'});
set(sgt, 'FontWeight', 'bold');
style_light(fs1, 13);
save_figure_png(fs1, dirV, 'eig_sweep_maxabs', ...
                struct('Width',15,'Height',7,'FontSize',15,'TitleSize',16));
close(fs1);

% ---- all six magnitudes ----------------------------------------------
%  Again two panels, for the same reason. On the left all six are flat
%  lines. On the right each one is divided by its own value at i = 90 deg,
%  which turns a flat line into the actual per-eigenvalue variation.
fs2 = figure('Visible', vis, 'Color', 'w', 'Name', 'eig sweep all');

subplot(1,2,1);
semilogy(INC_EIG, absev_J, '-', 'LineWidth', 1.8); hold on;
yline(1, 'k:', 'LineWidth', 1.4, 'HandleVisibility', 'off');
xline(c.i_crit_deg, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
grid on; xlim([-10 190]); xticks(0:45:180);
xlabel('inclination  [deg]'); ylabel('|\lambda|');
title('raw magnitudes, log scale');

subplot(1,2,2);
ref = absev_J(i90,:);
plot(INC_EIG, absev_J ./ ref, '-', 'LineWidth', 1.8); hold on;
xline(c.i_crit_deg, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(c.i_crit_retro, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
grid on; xlim([-10 190]); xticks(0:45:180);
xlabel('inclination  [deg]');
ylabel('|\lambda(i)| / |\lambda(90^\circ)|');
title('each one normalised to its polar value');

sgt = sgtitle({'All six eigenvalue magnitudes of \Phi(T) at \epsilon = 1', ...
       'two stay on the unit circle; the other four form reciprocal pairs'});
set(sgt, 'FontWeight', 'bold');
style_light(fs2, 13);
save_figure_png(fs2, dirV, 'eig_sweep_allmag', ...
                struct('Width',15,'Height',7,'FontSize',15,'TitleSize',16));
close(fs2);

% ---- secular drift: the plot that answers the guide's question -------
fs3 = figure('Visible', vis, 'Color', 'w', 'Name', 'drift sweep');
hw = plot(iv, [res.d_argp], 'o-', 'Color', [0.85 0.10 0.10], ...
          'LineWidth', 2.2, 'MarkerSize', 8);
hold on;
ho = plot(iv, [res.d_RAAN], 's-', 'Color', [0.00 0.35 0.85], ...
          'LineWidth', 2.2, 'MarkerSize', 8);
yline(0, 'k-', 'LineWidth', 1.0, 'HandleVisibility', 'off');
hc = xline(c.i_crit_deg, 'k--', 'LineWidth', 1.4);
h2 = xline(c.i_crit_retro, 'k--', 'LineWidth', 1.4);
set(h2, 'HandleVisibility', 'off');
grid on; xlim([0 180]); xticks(0:15:180);
xlabel('inclination  [deg]');
ylabel(sprintf('drift over %d revolutions  [deg]', N_REVS));
legend([hw ho hc], {'\omega  (apsides)', '\Omega  (node)', ...
       'critical inclination'}, 'Location', 'best');
title({sprintf('Secular drift over %d revolutions', N_REVS), ...
       '\omega crosses zero at the critical inclination; \Omega does not'});
style_light(fs3, 14);
save_figure_png(fs3, dirV, 'drift_sweep', struct('Width',12,'Height',7));
close(fs3);

% ---- closure ---------------------------------------------------------
fs4 = figure('Visible', vis, 'Color', 'w', 'Name', 'closure sweep');
hb = semilogy(iv, [res.clos_min], 'o-', 'Color', [0.20 0.55 0.20], ...
              'LineWidth', 2.2, 'MarkerSize', 8);
hold on;
he = semilogy(iv, [res.clos_end], 's--', 'Color', [0.45 0.10 0.60], ...
              'LineWidth', 1.8, 'MarkerSize', 8);
hc = xline(c.i_crit_deg, 'k--', 'LineWidth', 1.4);
grid on; xlim([0 180]); xticks(0:15:180);
xlabel('inclination  [deg]');
ylabel('distance from the starting point  [km]');
legend([hb he hc], {'best over all revolutions', 'after the last revolution', ...
       'critical inclination'}, 'Location', 'best');
title({'Closure test', ...
       'no inclination returns the satellite to its starting state'});
style_light(fs4, 14);
save_figure_png(fs4, dirV, 'closure_sweep', struct('Width',12,'Height',7));
close(fs4);

%% =====================================================================
%  6. DATA OUT
%  =====================================================================

Tb = table(iv.', [res.maxabs_K].', [res.maxabs_J].', [res.symp_J].', ...
           [res.d_argp].', [res.d_argp_pred].', ...
           [res.d_RAAN].', [res.d_RAAN_pred].', ...
           [res.clos_min].', [res.clos_min_rev].', [res.clos_end].', ...
           [res.equatorial].', ...
    'VariableNames', {'i_deg','maxabs_lambda_kepler','maxabs_lambda_j2', ...
                      'symplecticity_j2','d_argp_deg','d_argp_brouwer_deg', ...
                      'd_RAAN_deg','d_RAAN_brouwer_deg', ...
                      'closure_best_km','closure_best_rev','closure_final_km', ...
                      'equatorial_branch'});
writetable(Tb, fullfile(dirD, 'inclination_sweep_results.csv'));

Te = table(INC_EIG.', maxabs_K, maxabs_J, symp_K, symp_J, jres_K, rank_K, ...
    'VariableNames', {'i_deg','maxabs_kepler','maxabs_j2', ...
                      'symp_kepler','symp_j2','jordan_resid_kepler', ...
                      'rank_PhiminusI_kepler'});
writetable(Te, fullfile(dirD, 'eigenvalue_sweep_fine.csv'));

save(fullfile(dirD, 'inclination_sweep_results.mat'), ...
     'res', 'INC_EIG', 'maxabs_K', 'maxabs_J', 'absev_K', 'absev_J', ...
     'symp_K', 'symp_J', 'jres_K', 'rank_K', ...
     's', 'c', 'N_REVS', 'TOL_TRAJ', 'TOL_STM');

fprintf('=========================================================\n');
disp(Tb);
fprintf('=========================================================\n');
fprintf(' figures : %s\n', OUT_ROOT);
fprintf(' data    : %s\n', dirD);
fprintf('=========================================================\n');

if QUICK_TEST
    fprintf('\n*** That was QUICK_TEST. Check the folders, then set\n');
    fprintf('    QUICK_TEST = false at the top and run again. ***\n');
end

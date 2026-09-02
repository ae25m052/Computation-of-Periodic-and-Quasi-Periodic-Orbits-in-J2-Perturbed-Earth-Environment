%% KEPLERIAN vs J2 ORBIT-TO-ORBIT COMPARISON
%  Same initial state (from Keplerian elements) is propagated twice:
%       epsilon = 0  ->  pure two-body   ->  Kepler orbit
%       epsilon = 1  ->  J2 included     ->  perturbed orbit

clear; clc; close all;

%% ---------------------------------------------------------------
%  1. Constants and non-dimensionalization
%  -----------------------------------------------------------------
mu_dim = 398600.4418;      % Earth GM [km^3/s^2]
Re_dim = 6378.1363;        % Earth equatorial radius [km]
J2     = 1.08262668e-3;    % Earth J2

DU = Re_dim;                       % length unit [km]
TU = sqrt(DU^3/mu_dim);            % time unit [s]
VU = DU/TU;                        % velocity unit [km/s]

mu_nd = 1;
Re_nd = Re_dim/DU;                 


%% ===============================================================
%  2. RUN SETTINGS  --  *** CHANGE THE NUMBER OF ORBITS HERE ***
%  ===============================================================

N_ORBITS = 60;          % <<<<<<  NUMBER OF ORBITS TO PROPAGATE  >>>>>>

SAMPLES_PER_ORBIT = 400;   

THETA_MODE = 'trueanomaly';   

TIME_UNIT = 'auto';    

WRAP_THETA = true;     

%  ---------------------------------------------------------------
COL_KEP = [0.00 0.35 0.85];   
COL_J2  = [0.85 0.10 0.10];   
COL_DR  = [0.45 0.10 0.60];   
COL_DTH = [0.00 0.45 0.45];   

FS = 12;                     

%  ---------------------------------------------------------------
%  Initial osculating elements (the "initial state" of the sketch)
%  ---------------------------------------------------------------
a_dim = 7500;              % semi-major axis [km]
e0    = 0.05;              % eccentricity
i0    = deg2rad(45);       % inclination
RAAN0 = deg2rad(30);       % RAAN
argp0 = deg2rad(20);       % argument of perigee
nu0   = deg2rad(0);        % true anomaly at epoch


%% ---------------------------------------------------------------
%  3. Build the (single, shared) initial condition
%  -----------------------------------------------------------------
[r0_dim, v0_dim] = coe2rv(a_dim, e0, i0, RAAN0, argp0, nu0, mu_dim);
x0_nd = [r0_dim/DU; v0_dim/VU];

T_dim = 2*pi*sqrt(a_dim^3/mu_dim);   % nominal Keplerian period [s]
T_nd  = T_dim/TU;

fprintf('--- Run settings ---\n');
fprintf('N_ORBITS          : %d\n', N_ORBITS);
fprintf('Samples per orbit : %d  (total %d)\n', SAMPLES_PER_ORBIT, ...
        N_ORBITS*SAMPLES_PER_ORBIT);
fprintf('Seed orbit        : a=%.1f km, e=%.3f, i=%.1f deg\n', ...
        a_dim, e0, rad2deg(i0));
fprintf('Nominal period T  : %.3f s (%.3f min)\n\n', T_dim, T_dim/60);


%% ---------------------------------------------------------------
%  4. Propagate both cases over the SAME time span
%  -----------------------------------------------------------------
tspan = linspace(0, N_ORBITS*T_nd, N_ORBITS*SAMPLES_PER_ORBIT);
opts  = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% epsilon = 0 : pure Keplerian
[tK, XK] = ode89(@(t,x) twobody_eom(t, x, mu_nd, 0, J2, Re_nd), ...
                  tspan, x0_nd, opts);
XK = XK.';

% epsilon = 1 : J2 included
[tJ, XJ] = ode89(@(t,x) twobody_eom(t, x, mu_nd, 1, J2, Re_nd), ...
                  tspan, x0_nd, opts);
XJ = XJ.';

tK_sec = tK.' * TU;
tJ_sec = tJ.' * TU;
total_sec = max(tK_sec(end), tJ_sec(end));

unitChoice = lower(TIME_UNIT);
if strcmp(unitChoice, 'auto')
    if     total_sec <  3*3600,  unitChoice = 'min';
    elseif total_sec < 5*86400,  unitChoice = 'hr';
    else,                        unitChoice = 'day';
    end
end
switch unitChoice
    case 'sec', tScale = 1;      tLabel = 'time  [s]';
    case 'min', tScale = 60;     tLabel = 'time  [min]';
    case 'hr',  tScale = 3600;   tLabel = 'time  [hours]';
    case 'day', tScale = 86400;  tLabel = 'time  [days]';
    otherwise,  error('TIME_UNIT must be auto, sec, min, hr or day.');
end
tK_plot = tK_sec / tScale;
tJ_plot = tJ_sec / tScale;

fprintf('Time axis unit    : %s  (run length %.2f %s)\n\n', ...
        unitChoice, total_sec/tScale, unitChoice);


%% ---------------------------------------------------------------
%  5. Extract r(t), theta(t) and the revolution index
%  -----------------------------------------------------------------
[rK_nd, nuK, uK] = orbit_r_theta(XK, mu_nd);
[rJ_nd, nuJ, uJ] = orbit_r_theta(XJ, mu_nd);

rK = rK_nd * DU;      
rJ = rJ_nd * DU;

switch lower(THETA_MODE)
    case 'trueanomaly', thK = nuK;  thJ = nuJ;  thName = '\nu (true anomaly)';
    case 'arglat',      thK = uK;   thJ = uJ;   thName = 'u (arg. of latitude)';
    otherwise, error('THETA_MODE must be ''trueanomaly'' or ''arglat''.');
end

thK_un = unwrap(thK);
thJ_un = unwrap(thJ);

if WRAP_THETA
    thK_plot = rad2deg(mod(thK, 2*pi));
    thJ_plot = rad2deg(mod(thJ, 2*pi));
    thLabel  = sprintf('%s  [deg]', thName);
else
    thK_plot = rad2deg(thK_un);
    thJ_plot = rad2deg(thJ_un);
    thLabel  = sprintf('%s unwrapped  [deg]', thName);
end

uJ_un  = unwrap(uJ);
revJ   = floor((uJ_un - uJ_un(1)) / (2*pi)) + 1;   
uK_un  = unwrap(uK);
revK   = floor((uK_un - uK_un(1)) / (2*pi)) + 1;

nRevJ = max(revJ);
nRevK = max(revK);

fprintf('--- Revolutions actually completed in the same elapsed time ---\n');
fprintf('Keplerian : %d\n', nRevK);
fprintf('J2        : %d\n', nRevJ);
fprintf('(J2 can show one extra partial revolution: its nodal period is\n');
fprintf(' slightly shorter than the nominal Keplerian period.)\n\n');

fprintf('--- Radius ---\n');
fprintf('Kepler  r range : %.3f  to  %.3f km\n', min(rK), max(rK));
fprintf('J2      r range : %.3f  to  %.3f km\n', min(rJ), max(rJ));
fprintf('max |r_J2 - r_Kep| over the run : %.4f km\n\n', max(abs(rJ - rK)));


%% ---------------------------------------------------------------
%  6. FIGURE 1 - Keplerian case only
%  -----------------------------------------------------------------
figure('Name','Keplerian: r and theta vs time','Color','w');

subplot(2,1,1);
plot(tK_plot, rK, '-', 'Color', COL_KEP, 'LineWidth', 1.4); grid on;
ylabel('r  [km]'); set(gca,'GridAlpha',0.3);
title(sprintf('Keplerian (\\epsilon = 0):  radius over %d orbits', N_ORBITS));

subplot(2,1,2);
plot(tK_plot, thK_plot, '-', 'Color', COL_KEP, 'LineWidth', 1.4); grid on;
xlabel(tLabel); ylabel(thLabel); set(gca,'GridAlpha',0.3);
title('Keplerian (\epsilon = 0):  angle over time');

style_light(gcf, FS);


%% ---------------------------------------------------------------
%  7. FIGURE 2 - J2 case only
%  -----------------------------------------------------------------
figure('Name','J2: r and theta vs time','Color','w');

subplot(2,1,1);
plot(tJ_plot, rJ, '-', 'Color', COL_J2, 'LineWidth', 1.4); grid on;
ylabel('r  [km]'); set(gca,'GridAlpha',0.3);
title(sprintf('J_2-perturbed (\\epsilon = 1):  radius over %d orbits', N_ORBITS));

subplot(2,1,2);
plot(tJ_plot, thJ_plot, '-', 'Color', COL_J2, 'LineWidth', 1.4); grid on;
xlabel(tLabel); ylabel(thLabel); set(gca,'GridAlpha',0.3);
title('J_2-perturbed (\epsilon = 1):  angle over time');

style_light(gcf, FS);


%% ---------------------------------------------------------------
%  8. FIGURE 3 - Both cases on the SAME axes
%  -----------------------------------------------------------------
figure('Name','Kepler vs J2 overlaid','Color','w');

subplot(2,1,1);
plot(tK_plot, rK, '-', 'Color', COL_KEP, 'LineWidth', 2.2); hold on;
plot(tJ_plot, rJ, '-', 'Color', COL_J2, 'LineWidth', 1.2);
grid on; set(gca,'GridAlpha',0.3);
ylabel('r  [km]');
legend('Keplerian (\epsilon=0)', 'J_2 (\epsilon=1)', 'Location','best');
title(sprintf('Radius: Keplerian vs J_2  (%d orbits)', N_ORBITS));

subplot(2,1,2);
plot(tK_plot, thK_plot, '-', 'Color', COL_KEP, 'LineWidth', 2.2); hold on;
plot(tJ_plot, thJ_plot, '-', 'Color', COL_J2, 'LineWidth', 1.2);
grid on; set(gca,'GridAlpha',0.3);
xlabel(tLabel); ylabel(thLabel);
legend('Keplerian (\epsilon=0)', 'J_2 (\epsilon=1)', 'Location','best');
title('Angle: Keplerian vs J_2');

style_light(gcf, FS);


%% ---------------------------------------------------------------
%  9. FIGURE 4 - Orbit-to-orbit difference
%  -----------------------------------------------------------------
figure('Name','J2 minus Kepler','Color','w');

subplot(2,1,1);
plot(tK_plot, rJ - rK, '-', 'Color', COL_DR, 'LineWidth', 1.6); grid on;
ylabel('\Deltar  [km]'); set(gca,'GridAlpha',0.3);
title('Difference:  r_{J2} - r_{Kepler}');

subplot(2,1,2);
plot(tK_plot, rad2deg(thJ_un - thK_un), '-', 'Color', COL_DTH, 'LineWidth', 1.6); grid on;
xlabel(tLabel); ylabel('\Delta\theta  [deg]');
set(gca,'GridAlpha',0.3);
title('Difference:  \theta_{J2} - \theta_{Kepler}  (unwrapped)');

style_light(gcf, FS);


%% ---------------------------------------------------------------
% 10. FIGURE 5 - 3D J2 trajectory, ONE COLOUR PER REVOLUTION
%  -----------------------------------------------------------------
figure('Name','J2 trajectory coloured per orbit','Color','w');
plot_earth(Re_nd); hold on;

cmap = turbo(max(nRevJ,2));      
hRev = gobjects(1, nRevJ);       

for k = 1:nRevJ
    idx = find(revJ == k);
    if isempty(idx), continue; end
    if idx(end) < numel(revJ)
        idx = [idx, idx(end)+1]; 
    end
    hRev(k) = plot3(XJ(1,idx), XJ(2,idx), XJ(3,idx), '-', ...
                    'Color', cmap(k,:), 'LineWidth', 1.6);
end

grid on; set(gca,'GridAlpha',0.25); view(45,25); axis vis3d;
xlabel('x [DU]'); ylabel('y [DU]'); zlabel('z [DU]');
title(sprintf('J_2-perturbed trajectory: %d revolutions, one colour each', nRevJ));
lim = 1.15*max(vecnorm(XJ(1:3,:)));
xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);

if nRevJ <= 8
    valid = isgraphics(hRev);
    legend(hRev(valid), ...
           arrayfun(@(k) sprintf('Orbit %d',k), find(valid), 'UniformOutput', false), ...
           'Location','bestoutside');
else
    colormap(gca, cmap);
    cb = colorbar; cb.Label.String = 'orbit number';
    clim_safe(gca, [0.5, nRevJ+0.5]);
end

style_light(gcf, FS);


%% ---------------------------------------------------------------
% 11. FIGURE 6 - 3D Keplerian trajectory (reference, single colour)
%  -----------------------------------------------------------------
figure('Name','Keplerian trajectory','Color','w');
plot_earth(Re_nd); hold on;
plot3(XK(1,:), XK(2,:), XK(3,:), '-', 'Color', [0.20 0.85 0.95], 'LineWidth', 1.6);
grid on; set(gca,'GridAlpha',0.25); view(45,25); axis vis3d;
xlabel('x [DU]'); ylabel('y [DU]'); zlabel('z [DU]');
title(sprintf('Keplerian trajectory: %d orbits retrace one closed ellipse', N_ORBITS));
lim = 1.15*max(vecnorm(XK(1:3,:)));
xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);

style_light(gcf, FS);


%% ===============================================================
%  Local functions
%  ===============================================================
function [r, nu, u] = orbit_r_theta(X, mu)
%ORBIT_R_THETA Radius, true anomaly and argument of latitude along a run.
%   X  - 6xN state history (non-dimensional)
%   mu - gravitational parameter (matching units)
%   r  - 1xN radius
%   nu - 1xN true anomaly           [rad], angle measured from perigee
%   u  - 1xN argument of latitude   [rad], angle measured from the
%        ascending node. Well conditioned even as e -> 0, which is why the
%        revolution counter is built from u rather than from nu.

    N  = size(X,2);
    r  = zeros(1,N);  nu = zeros(1,N);  u = zeros(1,N);

    for k = 1:N
        rv = X(1:3,k);  vv = X(4:6,k);
        rn = norm(rv);
        r(k) = rn;

        h  = cross(rv,vv);   hn = norm(h);
        inc = acos(h(3)/hn);

        nvec = cross([0;0;1], h);  nn = norm(nvec);

        % --- argument of latitude ---
        cosu = dot(nvec, rv)/(nn*rn);
        sinu = rv(3)/(rn*sin(inc));
        u(k) = atan2(sinu, min(max(cosu,-1),1));

        % --- true anomaly ---
        evec = cross(vv,h)/mu - rv/rn;
        en   = norm(evec);
        cosnu = dot(evec, rv)/(en*rn);
        nuk   = acos(min(max(cosnu,-1),1));
        if dot(rv,vv) < 0, nuk = 2*pi - nuk; end
        nu(k) = nuk;
    end
end


function clim_safe(ax, lims)
    if exist('clim','file') || exist('clim','builtin')
        clim(ax, lims);
    else
        caxis(ax, lims); 
    end
end

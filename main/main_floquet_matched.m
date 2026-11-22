% Floquet-matched deputy on the circular chief orbits, h = 300 km, i = 1..89 deg.
% Deputy starts 200 m higher at the node. 

clear; clc; close all

c  = j2_constants();
mu = 1;  Re = 1;  J2 = c.J2;                   % canonical units
tod = @(W) rad2deg(W/c.TU*86400);              % rad/TU -> deg/day

alt  = 300;                                    % altitude at the node, km
inc  = 1:89;                                   % inclination, deg
dr   = 0.2/c.DU;                               % deputy 200 m higher
tend = 0.5*365.25*86400/c.TU;                  % half a year, TU
file = 'floquet_matched_results.xlsx';

r0 = (c.Re + alt)/c.DU;
n  = sqrt(mu/r0^3);
Vk = sqrt(mu/r0);

names = {'inc_deg','V_kmps','Wf_degday','T_min','laps','theta_radlap','cycle_days', ...
         'maxabs_lambda','dV_mps','di_deg','G_floq_km','G_nl_km', ...
         'late_floq_km','late_nl_km','diff_m'};
Res = NaN(numel(inc), numel(names));
Run = cell(numel(inc),1);                      % lap-by-lap data

fprintf(' i    theta [rad/lap]  cycle [d]  dV [m/s]     di [deg]      G floq/nl [km]\n');
for k = 1:numel(inc)
    i  = deg2rad(inc(k));
    Wb = -1.5*n*J2*(Re/r0)^2*cos(i);           % Brouwer node rate = first guess

    % chief, DP and cT
    [u,V,W,T,~,DP,cT] = shoot_circular(r0,i,0,Vk,Wb,mu,Re,J2,c.opts);
    X0c = [r0; 0; 0; u; V*cos(i) - W*r0; V*sin(i)];
    vc  = X0c(4:6);                            % chief velocity at the node
    N   = floor(tend/T);                       % laps in half a year

    % eigenvalues of DP, oscillation pair exp(+-i theta)
    lam = eig(DP);
    lam = lam(abs(lam) > 0.5);                 % drop the eigenvalue 0 (flow direction)
    [~,m] = max(abs(angle(lam)));
    th  = abs(angle(lam(m)));                  % perigee turn per lap, rad
    Ncy = round(2*pi/th);                      % laps in one perigee cycle

    row = [inc(k), V*c.VU, tod(W), T*c.TU/60, N, th, Ncy*T*c.TU/86400, ...
           max(abs(lam)), NaN(1,7)];

    % option A needs at least one whole perigee cycle in half a year
    if Ncy <= N
        % deputy start for a speed change dV and tilt change di
        dX0 = @(dV,di) [r0+dr; 0; 0; u; (V+dV)*cos(i+di) - W*(r0+dr); (V+dV)*sin(i+di)] - X0c;

        % Newton on dV, di: zero average dy growth and zero average dT
        Nm = floor(N/Ncy)*Ncy;                 % whole cycles
        q  = [0; 0];
        h  = 1e-9;
        for s = 1:6
            F  = secular(dX0(q(1),q(2)), DP, cT, Nm);
            Jq = [secular(dX0(q(1)+h,q(2)), DP, cT, Nm) - F, ...
                  secular(dX0(q(1),q(2)+h), DP, cT, Nm) - F]/h;
            q  = q - Jq\F;
        end
        d = dX0(q(1),q(2));

        % Floquet prediction and nonlinear run
        [dXf,dTf] = floquet_laps(d, DP, cT, N);
        [dXn,dTn] = nonlinear_laps(X0c, d, W, T, N, c);
        sf = same_time(dXf, dTf, vc);
        sn = same_time(dXn, dTn, vc);
        Run{k} = struct('dXf',dXf,'dTf',dTf,'dXn',dXn,'dTn',dTn,'sf',sf,'sn',sn,'T',T, ...
                        'X0c',X0c,'W',W,'dX0',d,'N',N);

        row(9:15) = [q(1)*c.VU*1000, rad2deg(q(2)), ...
                     (max(sf) - sf(1))*c.DU, (max(sn) - sn(1))*c.DU, ...
                     V*sum(dTf)*c.DU, V*sum(dTn)*c.DU, ...
                     max(vecnorm(dXf(1:3,:) - dXn(1:3,:)))*c.DU*1000];
    end
    Res(k,:) = row;

    fprintf('%2d   %.6f       %7.1f    %+.6f    %+.3e    %8.3f %8.3f\n', ...
            inc(k), th, row(7), row(9), row(10), row(11), row(12));
end

%% Excel file
if isfile(file), delete(file); end
writetable(array2table(Res, 'VariableNames', names), file, 'Sheet', 'Matched');
fprintf('\nResults written to %s\n', file);

%% Which chief is good (smallest growth, nonlinear)
[~,kb] = min(Res(:,12));
[~,kw] = max(Res(:,12));
fprintf('Best chief:  i = %d deg, G = %.3f km\n', inc(kb), Res(kb,12));
fprintf('Worst chief: i = %d deg, G = %.3f km\n', inc(kw), Res(kw,12));
fprintf('Option A not possible (cycle longer than half a year): i = %s deg\n', ...
        mat2str(inc(isnan(Res(:,9)))));

figure; hold on
plot(inc, Res(:,11), 'b-', 'DisplayName', 'Floquet DP^k')
plot(inc, Res(:,12), 'ro', 'DisplayName', 'nonlinear')
grid on; xlabel('i  [deg]'); ylabel('G  [km]')
title('Floquet matched \delta X_0: growth in half a year'); legend('Location','best')
fig_format(gcf)

%% Floquet numbers against inclination
figure
subplot(2,2,1)
plot(inc, Res(:,6), 'b.-'); grid on
xlabel('i  [deg]'); ylabel('\theta  [rad/lap]'); title('angle of e^{\pmi\theta}')
subplot(2,2,2)
semilogy(inc, Res(:,7), 'b.-'); hold on
yline(tend*c.TU/86400, 'r--', 'half a year'); grid on
xlabel('i  [deg]'); ylabel('perigee cycle  [days]'); title('2\pi/\theta laps')
subplot(2,2,3)
plot(inc, Res(:,9)*1000, 'bo'); grid on
xlabel('i  [deg]'); ylabel('\delta V  [mm/s]'); title('deputy speed change')
subplot(2,2,4)
plot(inc, Res(:,10), 'bo'); grid on
xlabel('i  [deg]'); ylabel('\delta i  [deg]'); title('deputy tilt change')
fig_format(gcf, 30)

%% Lap-by-lap offsets: best, worst and 45 deg
show = unique([kb, kw, find(inc == 45 & ~isnan(Res(:,9)).')]);
lb = {'x','y','z'};
for k = show
    R   = Run{k};
    lap = 0:R.N;
    dd  = lap*R.T*c.TU/86400;                  % days
    figure('Name', sprintf('i = %d deg', inc(k)))
    for j = 1:3
        subplot(2,2,j); hold on
        plot(lap, R.dXf(j,:)*c.DU*1000, 'b-', 'DisplayName', 'Floquet DP^k')
        plot(lap, R.dXn(j,:)*c.DU*1000, 'r--', 'DisplayName', 'nonlinear')
        grid on; xlabel('lap'); ylabel(['\delta ' lb{j} '  [m]'])
        title(['\delta ' lb{j} ' at the node'])
    end
    legend('Location','best')
    subplot(2,2,4); hold on
    plot(dd, R.sf*c.DU, 'b-', 'DisplayName', 'Floquet DP^k')
    plot(dd, R.sn*c.DU, 'r--', 'DisplayName', 'nonlinear')
    grid on; xlabel('t  [days]'); ylabel('separation  [km]')
    title('same-time separation')
    sgtitle(sprintf('h = %d km, i = %d deg, Floquet matched \\delta X_0, frame turning at W_f', alt, inc(k)))
    fig_format(gcf, 30)

    % relative orbit over half a year, frame turning at W_f
    t  = linspace(0, R.N*R.T, 20*R.N + 1);
    Y0 = [R.X0c; R.X0c + R.dX0];
    [t,Y] = ode89(@(t,Y) [rf_eom(t,Y(1:6),mu,Re,J2,R.W); rf_eom(t,Y(7:12),mu,Re,J2,R.W)], ...
                  t, Y0, c.opts);
    dP = (Y(:,7:9) - Y(:,1:3)).'*c.DU*1000;    % m
    figure; hold on
    plot_orbit_gradient(dP, t*c.TU/86400, 'turbo', 1);
    cb = colorbar; cb.Label.String = 't  [days]';
    grid on; axis equal tight; view(45,25)
    xlabel('\delta x  [m]'); ylabel('\delta y  [m]'); zlabel('\delta z  [m]')
    title(sprintf('Deputy minus chief, i = %d deg, half a year', inc(k)))
    fig_format(gcf)
end

%% Best chief: Earth view (frame turning at W_f) and LVLH zoom, saved as a GIF
R  = Run{kb};
t  = linspace(0, R.N*R.T, 20*R.N + 1);         % 20 points per lap
Y0 = [R.X0c; R.X0c + R.dX0];
[t,Y] = ode89(@(t,Y) [rf_eom(t,Y(1:6),mu,Re,J2,R.W); rf_eom(t,Y(7:12),mu,Re,J2,R.W)], ...
              t, Y0, c.opts);
rc = Y(:,1:3).';  rd = Y(:,7:9).';             % chief, deputy positions
vI = Y(:,4:6).' + cross(repmat([0; 0; R.W],1,numel(t)), rc);   % chief inertial velocity
L  = lvlh(rc, vI, rd - rc)*c.DU*1000;          % deputy in chief LVLH, m
P  = [L(2,:); L(1,:); L(3,:)];                 % plot order: along-track, radial, cross-track
td = t*c.TU/86400;                             % days
[~,Xc] = ode89(@(t,X) rf_eom(t,X,mu,Re,J2,R.W), linspace(0,R.T,400), R.X0c, c.opts);   % one chief lap

fig = figure('Name', sprintf('Best chief, i = %d deg', inc(kb)));
ax1 = subplot(1,2,1); hold on
he = plot_earth(1);
set(he, 'HandleVisibility', 'off')
h1 = plot3(Xc(:,1), Xc(:,2), Xc(:,3), 'k', 'LineWidth', 1.5);
hc = plot3(rc(1,1), rc(2,1), rc(3,1), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 9);
hd = plot3(rd(1,1), rd(2,1), rd(3,1), 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
lim = 1.3;  xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim])
grid on; axis vis3d; view(45,25)
xlabel('x [DU]'); ylabel('y [DU]'); zlabel('z [DU]')
title('Earth view, frame turning at W_f')
legend([h1 hc hd], {'chief orbit','chief','deputy (0.2-0.6 km away)'}, 'Location', 'southoutside')

ax2 = subplot(1,2,2); hold on
hg = plot_orbit_gradient(P(:,1), td(1), 'turbo', 1);
clim(ax2, [0 td(end)])
cb = colorbar; cb.Label.String = 't  [days]';
plot3(0, 0, 0, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 12)
hl = plot3(P(1,1), P(2,1), P(3,1), 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
xlim(1.05*[min(P(1,:)) max(P(1,:))]); ylim(1.05*[min(P(2,:)) max(P(2,:))])
zlim(max(1.5*max(abs(P(3,:))), 1)*[-1 1])
grid on; view(-30,25)
xlabel('along-track [m]'); ylabel('radial [m]'); zlabel('cross-track [m]')
title('Deputy seen from the chief (LVLH)')
ht = sgtitle(sprintf('Best chief i = %d deg, Floquet matched deputy, day 0', inc(kb)));
fig_format(fig, 34)

% GIF: each frame jumps 15 laps + 1/20 lap, so the satellites creep round the orbit,
% and the LVLH trail is drawn up to the current time
gif = sprintf('best_chief_i%d.gif', inc(kb));
for j = 1:301:numel(t)
    set(hc, 'XData', rc(1,j), 'YData', rc(2,j), 'ZData', rc(3,j))
    set(hd, 'XData', rd(1,j), 'YData', rd(2,j), 'ZData', rd(3,j))
    delete(hg)
    set(fig, 'CurrentAxes', ax2)
    hg = plot_orbit_gradient(P(:,1:j), td(1:j), 'turbo', 1);
    clim(ax2, [0 td(end)])
    set(hl, 'XData', P(1,j), 'YData', P(2,j), 'ZData', P(3,j))
    uistack(hl, 'top')
    ht.String = sprintf('Best chief i = %d deg, Floquet matched deputy, day %.0f, distance %.0f m', ...
                        inc(kb), td(j), norm(L(:,j)));
    drawnow
    [im,map] = rgb2ind(frame2im(getframe(fig)), 256);
    if j == 1
        imwrite(im, map, gif, 'gif', 'LoopCount', inf, 'DelayTime', 0.05);
    else
        imwrite(im, map, gif, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end
fprintf('GIF written to %s\n', gif);


function [dX,dT] = floquet_laps(d,DP,cT,N)
% dX(:,k+1) = DP^k d (offset at the k-th node), dT(k) = lap-time change of lap k
dX = zeros(6,N+1);  dT = zeros(1,N);
dX(:,1) = d;
for k = 1:N
    dT(k)     = cT*dX(:,k);
    dX(:,k+1) = DP*dX(:,k);
end
end


function [dX,dT] = nonlinear_laps(X0c,d,W,T,N,c)
% Deputy from X0c + d with rf_eom, offsets at its first N upward nodes.
mu = 1;  Re = 1;  J2 = c.J2;
o = odeset(c.opts,'Events',@event_node_all);
[~,~,te,Xe] = ode89(@(t,X) rf_eom(t,X,mu,Re,J2,W), [0 (N + 0.5)*T], X0c + d, o);
m  = te > 1e-6*T;
te = te(m);  Xe = Xe(m,:);
te = te(1:N);  Xe = Xe(1:N,:);
dX = [d, Xe.' - X0c];
dT = diff([0; te]).' - T;
end


function F = secular(d,DP,cT,N)
% Average y growth per lap and average lap-time change over N laps.
[dX,dT] = floquet_laps(d,DP,cT,N);
p = polyfit(0:N, dX(2,:), 1);
F = [p(1); mean(dT)];
end


function s = same_time(dX,dT,vc)
% Separation at the same time as the chief's k-th node.
s = vecnorm(dX(1:3,:) - vc*[0 cumsum(dT)]);
end


function L = lvlh(r, v, d)
% Components of d along radial, along-track and cross-track of the orbit (r, v).
er = r./vecnorm(r);
eh = cross(r, v);  eh = eh./vecnorm(eh);
et = cross(eh, er);
L  = [sum(d.*er); sum(d.*et); sum(d.*eh)];
end


function [val,stop,dir] = event_node_all(~,X)
val  = X(3);     % z = 0
stop = 0;
dir  = 1;        % going up only
end

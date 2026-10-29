% J2 orbits made periodic in a frame spinning about z at rate Wf.
% Part 1: near-circular orbits at fixed inclinations.
% Part 2: eccentric orbit, inclination found by Newton.

clear; clc; close all

c = j2_constants();

% normalise first
mu = 1;
Re = 1;
J2 = c.J2;

hp = 300;                         % km
rp = (c.Re + hp)/c.DU;            % DU

%% Part 1: near-circular, fixed i
inc = [1 15 30 45 60 75 90 105 120 135 150 165 179];   % deg

N  = numel(inc);
WF = zeros(N,1);  WB = zeros(N,1);

fprintf('Near-circular orbits, hp = %.1f km\n',hp);
fprintf('  i [deg]  it   vp [km/s]       Wf [deg/day]  Brouwer [deg/day]  T [min]    closure [km]\n');

n = sqrt(mu/rp^3);

for k = 1:N
    i = deg2rad(inc(k));
    [vp,~,Wf,T,hist] = shoot_inclined(rp,[],i,mu,Re,J2,c.opts);

    X0 = [rp*cos(i); 0; rp*sin(i); 0; vp - Wf*rp*cos(i); 0];
    [~,X] = ode89(@(t,X) rf_eom(t,X,mu,Re,J2,Wf),[0 T],X0,c.opts);
    err = norm(X(end,1:3) - X0(1:3).');

    WF(k) = rad2deg(Wf/c.TU*86400);
    WB(k) = rad2deg(-1.5*n*J2*(Re/rp)^2*cos(i)/c.TU*86400);

    fprintf('  %6.1f   %d   %.10f   %+9.4f     %+9.4f        %.4f   %.2e\n', ...
            inc(k),size(hist,1)-1,vp*c.VU,WF(k),WB(k),T*c.TU/60,err*c.DU);
end

ic = 0:1:180;
WBc = rad2deg(-1.5*n*J2*(Re/rp)^2*cosd(ic)/c.TU*86400);

figure; hold on; grid on; box on
plot(ic,WBc,'r-','LineWidth',1.5)
plot(inc,WF,'bo','MarkerFaceColor','b','MarkerSize',6)
yline(0,'k:');
xlim([0 180]); xticks(0:30:180)
xlabel('inclination  [deg]'); ylabel('\Omega_f  [deg/day]')
title(sprintf('Frame rate vs inclination, near-circular, h_p = %.0f km',hp))
legend('Brouwer nodal rate','shooting','Location','northwest')
fig_format(gcf)

%% Part 2: eccentric, i free
e  = 0.65;
ra = rp*(1+e)/(1-e);

[vp,i,Wf,T,hist] = shoot_inclined(rp,ra,[],mu,Re,J2,c.opts);

a  = rp/(1-e);
p  = a*(1-e^2);
n  = sqrt(mu/a^3);
WB = -1.5*n*J2*(Re/p)^2*cos(i);

fprintf('\nEccentric orbit, hp = %.1f km, e = %.2f, %d Newton steps\n',hp,e,size(hist,1)-1);
fprintf('  i  = %.6f deg\n',rad2deg(i));
fprintf('  vp = %.10f km/s  (%+.4f m/s from Kepler)\n',vp*c.VU,(vp-sqrt(mu*(1+e)/rp))*c.VU*1000);
fprintf('  Wf = %+.4f deg/day   Brouwer = %+.4f deg/day\n', ...
        rad2deg(Wf/c.TU*86400),rad2deg(WB/c.TU*86400));
fprintf('  T  = %.10f TU = %.4f min = %.6f h\n',T,T*c.TU/60,T*c.TU/3600);

X0 = [rp*cos(i); 0; rp*sin(i); 0; vp - Wf*rp*cos(i); 0];
[~,X] = ode89(@(t,X) rf_eom(t,X,mu,Re,J2,Wf),[0 T],X0,c.opts);
err = norm(X(end,1:3) - X0(1:3).');
fprintf('  closure after one period, rotating frame = %.2e km\n',err*c.DU);

% 30 periods: rotating frame and ECI
nrev = 30;
[tr,Xr] = ode89(@(t,X) rf_eom(t,X,mu,Re,J2,Wf),[0 nrev*T],X0,c.opts);
X0i = [rp*cos(i); 0; rp*sin(i); 0; vp; 0];
[ti,Xi] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,1),[0 nrev*T],X0i,c.opts);

plot_orbit(tr,Xr,T,sprintf('Rotating frame (\\Omega_f = %+.4f deg/day), i = %.2f deg, 30 periods', ...
           rad2deg(Wf/c.TU*86400),rad2deg(i)));
fig_format(gcf)
plot_orbit(ti,Xi,T,sprintf('ECI, i = %.2f deg, 30 periods (node drifts)',rad2deg(i)));
fig_format(gcf)

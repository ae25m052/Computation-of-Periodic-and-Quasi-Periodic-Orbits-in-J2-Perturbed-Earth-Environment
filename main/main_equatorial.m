% Equatorial case, i = 0. The perigee speed is corrected by a Newton step
% so that the J2 orbit reaches the intended apogee radius.

clear; clc; close all

c = j2_constants();

% normalise first
mu = 1;
Re = 1;
J2 = c.J2;

hp = 300;                         % km
e  = 0.65;
rp = (c.Re + hp)/c.DU;            % DU
ra_target = rp*(1+e)/(1-e);

fprintf('rp = %.10f DU = %.4f km (altitude %.1f km)\n',rp,rp*c.DU,hp);
fprintf('target ra = %.10f DU = %.4f km\n\n',ra_target,ra_target*c.DU);

[vp,ra,Trad,hist] = shoot_equatorial(rp,ra_target,mu,Re,J2,c.opts);

fprintf('Newton iterations\n');
fprintf('  k        vp [DU/TU]          ra error [DU]       ra error [km]\n');
for k = 1:size(hist,1)
    fprintf('  %d   %.13f   %+14.6e   %+14.6e\n',hist(k,1),hist(k,2),hist(k,3),hist(k,3)*c.DU);
end

vp_kepler = sqrt(mu*(1+e)/rp);

fprintf('\nKepler perigee speed    = %.13f DU/TU = %.9f km/s\n',vp_kepler,vp_kepler*c.VU);
fprintf('corrected perigee speed = %.13f DU/TU = %.9f km/s\n',vp,vp*c.VU);
fprintf('difference              = %.6f m/s\n',(vp-vp_kepler)*c.VU*1000);
fprintf('\nradial period = %.10f TU = %.6f min\n',Trad,Trad*c.TU/60);
fprintf('Kepler period = %.10f TU = %.6f min\n',2*pi*sqrt((rp/(1-e))^3/mu), ...
        2*pi*sqrt((rp/(1-e))^3/mu)*c.TU/60);

% propagate the corrected orbit
X0 = [rp;0;0; 0;vp;0];
[~,X] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,1),[0 Trad],X0,c.opts);

figure; hold on; grid on; axis equal
th = linspace(0,2*pi,200);
fill(Re*cos(th),Re*sin(th),[0.25 0.45 0.75],'EdgeColor','none')
plot(X(:,1),X(:,2),'r','LineWidth',1.2)
xlabel('x [DU]'); ylabel('y [DU]')
title('Equatorial J2 orbit after correction')

% closure over one radial period
err = norm(X(end,1:3)-X0(1:3).');
fprintf('\nposition error after one radial period = %.6e DU = %.4f km\n',err,err*c.DU);

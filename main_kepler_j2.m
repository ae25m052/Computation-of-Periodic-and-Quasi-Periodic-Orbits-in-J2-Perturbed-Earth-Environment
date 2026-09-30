% Keplerian orbit and J2-perturbed orbit from the same initial state.

clear; clc; close all

c = j2_constants();

% normalise first
mu = 1;
Re = 1;
J2 = c.J2;

a    = 7500/c.DU;                 % DU
e    = 0.05;
i    = deg2rad(45);
RAAN = deg2rad(30);
argp = deg2rad(20);
nu   = 0;

[r0,v0] = kep2rv(a,e,i,RAAN,argp,nu,mu);
X0 = [r0; v0];

T = 2*pi*sqrt(a^3/mu);            % TU

fprintf('a = %.6f DU = %.1f km, e = %.2f, i = %.1f deg\n',a,a*c.DU,e,rad2deg(i));
fprintf('T = %.8f TU = %.4f s = %.4f min\n\n',T,T*c.TU,T*c.TU/60);

% one period, Kepler and J2
[tk,Xk] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,0),[0 T],X0,c.opts);
[tj,Xj] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,1),[0 T],X0,c.opts);

ek = norm(Xk(end,1:3)-X0(1:3).');
ej = norm(Xj(end,1:3)-X0(1:3).');
fprintf('closure error after one period\n');
fprintf('  Kepler : %.3e DU = %.3e km\n',ek,ek*c.DU);
fprintf('  J2     : %.6e DU = %.3f km\n\n',ej,ej*c.DU);

plot_orbit(tk,Xk,T,'Keplerian orbit');
plot_orbit(tj,Xj,T,'J2-perturbed orbit');

% 30 periods to show the drift
[t30,X30] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,1),[0 30*T],X0,c.opts);
plot_orbit(t30,X30,T,'J2 orbit, 30 periods');

% element history
% n = numel(t30);
% RA = zeros(n,1); AP = zeros(n,1);
% for k = 1:n
%     [~,~,~,RA(k),AP(k)] = rv2kep(X30(k,1:3).',X30(k,4:6).',mu);
% end
% RA = unwrap(RA); AP = unwrap(AP);
% 
% figure
% subplot(2,1,1)
% plot(t30/T,rad2deg(RA-RA(1))); grid on
% ylabel('\Delta\Omega [deg]')
% subplot(2,1,2)
% plot(t30/T,rad2deg(AP-AP(1))); grid on
% xlabel('orbits'); ylabel('\Delta\omega [deg]')
% 
% fprintf('over 30 periods\n');
% fprintf('  RAAN drift = %.4f deg\n',rad2deg(RA(end)-RA(1)));
% fprintf('  argp drift = %.4f deg\n',rad2deg(AP(end)-AP(1)));

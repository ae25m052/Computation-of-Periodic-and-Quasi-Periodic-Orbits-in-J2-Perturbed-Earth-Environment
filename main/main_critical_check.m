% Critical inclination check. 

clear; clc; close all

c = j2_constants();

% normalise first
mu = 1;
Re = 1;
J2 = c.J2;

hp = 300;                         % km
e  = 0.65;
rp = (c.Re + hp)/c.DU;            % DU
a  = rp/(1-e);
p  = a*(1-e^2);
n  = sqrt(mu/a^3);
T  = 2*pi*sqrt(a^3/mu);           % TU

nrev = 100;

% wdot = 0 when 5 cos^2 i - 1 = 0
i_eq = acos(sqrt(1/5));

[r0,v0] = kep2rv(a,e,i_eq,0,0,0,mu);

% output at every T, events at every perigee, half a period extra so the
% last perigee is not missed
opts  = odeset(c.opts,'Events',@event_perigee);
tspan = [(0:nrev)*T, (nrev+0.5)*T];
[~,X,te,Xe] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,1),tspan,[r0;v0],opts);

% drop an event at t = 0 if the solver reports it, keep nrev perigees
Xe = Xe(te > 1e-6,:);
te = te(te > 1e-6);
Xp = [[r0;v0].'; Xe(1:nrev,:)];

iT = zeros(nrev+1,1);  wT = zeros(nrev+1,1);
iP = zeros(nrev+1,1);  wP = zeros(nrev+1,1);
for k = 1:nrev+1
    [~,~,iT(k)] = rv2kep(X(k,1:3).',X(k,4:6).',mu);
    [~,~,iP(k)] = rv2kep(Xp(k,1:3).',Xp(k,4:6).',mu);
    wT(k) = 0.75*n*J2*(Re/p)^2*(5*cos(iT(k))^2 - 1);
    wP(k) = 0.75*n*J2*(Re/p)^2*(5*cos(iP(k))^2 - 1);
end
wT = rad2deg(wT*T);               % deg/rev
wP = rad2deg(wP*T);

Wd = rad2deg(-1.5*n*J2*(Re/p)^2*cos(i_eq)*T);

fprintf('hp = %.1f km, e = %.2f, i from equation = %.6f deg\n',hp,e,rad2deg(i_eq));
fprintf('node rate at this i = %.4f deg/rev (for comparison)\n',Wd);
fprintf('Kepler period T = %.4f min, perigee to perigee = %.4f min\n\n', ...
        T*c.TU/60,(te(2)-te(1))*c.TU/60);
fprintf('                       every Kepler period T                  at every perigee\n');
fprintf('  rev   i_eq [deg]   i [deg]     i - i_eq    wdot [deg/rev]   i [deg]     i - i_eq    wdot [deg/rev]\n');
for k = 1:nrev+1
    fprintf('  %3d   %.6f    %.6f   %+.6f   %+.3e       %.6f   %+.6f   %+.3e\n', ...
            k-1,rad2deg(i_eq),rad2deg(iT(k)),rad2deg(iT(k)-i_eq),wT(k), ...
            rad2deg(iP(k)),rad2deg(iP(k)-i_eq),wP(k));
end


function [value,isterminal,direction] = event_perigee(~,X)
% Perigee: r.v goes from negative to positive.

value      = X(1:3).'*X(4:6);
isterminal = 0;
direction  = 1;

end

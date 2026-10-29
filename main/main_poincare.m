% Periodic J2 orbits as fixed points of the Poincare map on the equator,
% then the same orbit shown in three frames:
%   ECI (not turning), ECEF (turns once per sidereal day),
%   and the frame turning at the node rate Wf found by the shooting.

clear; clc; close all

c  = j2_constants();
mu = 1;  Re = 1;  J2 = c.J2;               

WE  = 2*pi/86164.0905*c.TU;             
tod = @(W) rad2deg(W/c.TU*86400);         

%% Orbit 1: near-circular, hp = 300 km, i = 45 deg  (r0, i fixed)
r0 = (c.Re + 300)/c.DU;
i1 = deg2rad(45);
n1 = sqrt(mu/r0^3);
Wb1 = -1.5*n1*J2*(Re/r0)^2*cos(i1);         

[i1,u1,V1,W1,T1,h1] = shoot_poincare('circ',r0,i1,0,sqrt(mu/r0),Wb1,mu,Re,J2,c.opts);

fprintf('Orbit 1: near-circular, hp = 300 km, i = 45 deg\n');
fprintf('  step   u [m/s]      V [km/s]        Wf [deg/day]   |F| [DU]\n');
for k = 1:size(h1,1)
    fprintf('  %d     %+.3e   %.10f   %+.6f     %.1e\n', ...
        h1(k,1), h1(k,3)*c.VU*1000, h1(k,4)*c.VU, tod(h1(k,5)), h1(k,6));
end
fprintf('  T  = %.6f min (node to node)\n', T1*c.TU/60);
fprintf('  Wf = %+.6f deg/day,  Brouwer = %+.6f deg/day\n\n', tod(W1), tod(Wb1));

%% Orbit 2: eccentric, hp = 300 km, e = 0.65, perigee at the node  (r0, u, V fixed)
e  = 0.65;
a2 = r0/(1 - e);
p2 = a2*(1 - e^2);
n2 = sqrt(mu/a2^3);
V2 = sqrt(mu*(1 + e)/r0);                     
i2 = acos(sqrt(1/5)) + deg2rad(0.3);          
Wb2 = -1.5*n2*J2*(Re/p2)^2*cos(i2);

[i2,u2,V2,W2,T2,h2] = shoot_poincare('ecc',r0,i2,0,V2,Wb2,mu,Re,J2,c.opts);

fprintf('Orbit 2: hp = 300 km, e = %.2f, perigee at the ascending node\n', e);
fprintf('  step   i [deg]         Wf [deg/day]   |F| [DU]\n');
for k = 1:size(h2,1)
    fprintf('  %d     %.8f    %+.6f     %.1e\n', h2(k,1), rad2deg(h2(k,2)), tod(h2(k,5)), h2(k,6));
end
fprintf('  T  = %.6f h (node to node)\n', T2*c.TU/3600);
fprintf('  Wf = %+.6f deg/day,  Brouwer = %+.6f deg/day\n\n', tod(W2), tod(-1.5*n2*J2*(Re/p2)^2*cos(i2)));

% Newton convergence
figure
semilogy(h1(:,1),h1(:,6),'bo-'); hold on
semilogy(h2(:,1),h2(:,6),'rs-'); grid on
xlabel('Newton step'); ylabel('miss |F|  [DU]')
title('Differential correction on the Poincare map')
legend('orbit 1 (near-circular)','orbit 2 (e = 0.65)')
xticks(0:max([h1(:,1); h2(:,1)]))
fig_format(gcf)

%% Stability: monodromy matrix after one lap
X01 = [r0; 0; 0; u1; V1*cos(i1) - W1*r0; V1*sin(i1)];
X02 = [r0; 0; 0; u2; V2*cos(i2) - W2*r0; V2*sin(i2)];
orb = {X01, W1, T1, 'Orbit 1'; X02, W2, T2, 'Orbit 2'};
for k = 1:2
    Y0 = [orb{k,1}; reshape(eye(6),36,1); zeros(6,1)];
    [~,Y] = ode89(@(t,Y) rf_stm_eom(t,Y,mu,Re,J2,orb{k,2}),[0 orb{k,3}],Y0,c.opts);
    lam = eig(reshape(Y(end,7:42),6,6));
    fprintf('%s monodromy eigenvalues (|lambda|):\n', orb{k,4});
    fprintf('  %+.8f %+.8fi   (%.8f)\n', [real(lam) imag(lam) abs(lam)].');
    fprintf('  closure after one lap, all 6 numbers: %.1e\n\n', norm(Y(end,1:6).' - orb{k,1}));
end

%% Same orbits in three frames
fr = {'ECI', 'ECEF (sidereal day)', 'Node-rate frame'};
show_frames(X01, W1, T1, 16, 'Orbit 1', fr, [0 WE W1], c);
show_frames(X02, W2, T2, 6,  'Orbit 2', fr, [0 WE W2], c);

%% Poincare maps: r and rdot at every ascending node
% Converged orbit plus starts pushed outward by a few m/s.
kick = [0 5 10 20 40]/1000/c.VU;               % m/s -> DU/TU
poincare_map(X01, W1, T1, 900, kick, 'Orbit 1', c);
poincare_map(X02, W2, T2, 300, kick, 'Orbit 2', c);   


function show_frames(X0, W, T, nlap, name, fr, Wk, c)
% Same physical orbit in frames turning at Wk(1), Wk(2), Wk(3).
% All frames line up at t = 0, so the start position is the same.

mu = 1;  Re = 1;  J2 = c.J2;
r  = X0(1:3);
vI = X0(4:6) + cross([0; 0; W], r);         

o = odeset(c.opts,'Events',@event_node_all);
fprintf('%s, %d laps:\n', name, nlap);
fprintf('  frame                    node shift per lap [deg]   start-to-last-node distance [km]\n');

figure('Name',[name ' nodes'])
for k = 1:3
    X0k = [r; vI - cross([0; 0; Wk(k)], r)];  
    t   = linspace(0, nlap*T, 300*nlap);
    [t,X,te,Xe] = ode89(@(t,X) rf_eom(t,X,mu,Re,J2,Wk(k)), t, X0k, o);
    Xe = Xe(te > 1e-6*T,:);                    

    lon = rad2deg(unwrap(atan2([0; Xe(:,2)], [r(1); Xe(:,1)])));
    fprintf('  %-24s %+12.4f              %12.3f\n', fr{k}, mean(diff(lon)), ...
            norm(Xe(end,1:3).' - r)*c.DU);

    % orbit
    plot_orbit(t, X, T, sprintf('%s in %s, %d laps', name, fr{k}, nlap));
    fig_format(gcf)

    % ascending nodes (Poincare section points)
    figure(findobj('Name',[name ' nodes']))
    subplot(1,3,k)
    th = linspace(0,2*pi,200);
    plot(c.Re*cos(th), c.Re*sin(th), 'k-'); hold on
    plot(Xe(:,1)*c.DU, Xe(:,2)*c.DU, 'ro', 'MarkerFaceColor','r', 'MarkerSize',4)
    plot(r(1)*c.DU, r(2)*c.DU, 'gp', 'MarkerFaceColor','g', 'MarkerSize',12)
    axis equal; grid on
    lim = 1.15*r(1)*c.DU;  xlim([-lim lim]); ylim([-lim lim])
    xlabel('x  [km]'); ylabel('y  [km]')
    title(fr{k})
end
sgtitle(sprintf('%s: ascending nodes over %d laps (star = start)', name, nlap))
fig_format(gcf, 30)
fprintf('\n');

end


function poincare_map(X0, W, T, nlap, kick, name, c)
% Crossings of z = 0 going up, for nlap laps, plotted as
% (r - r0, rdot). r and rdot are the same in every frame.

mu = 1;  Re = 1;  J2 = c.J2;
r0 = norm(X0(1:3));
o  = odeset(c.opts,'Events',@event_node_all);

figure; hold on; grid on; box on
for k = 1:numel(kick)
    X0k = X0 + [0; 0; 0; kick(k); 0; 0];      
    [~,~,te,Xe] = ode89(@(t,X) rf_eom(t,X,mu,Re,J2,W), [0 nlap*T], X0k, o);
    Xe = Xe(te > 1e-6*T,:);
    r  = sqrt(sum(Xe(:,1:3).^2,2));
    rd = sum(Xe(:,1:3).*Xe(:,4:6),2)./r;      
    if k == 1
        plot((r - r0)*c.DU, rd*c.VU*1000, 'ko', 'MarkerFaceColor','g', 'MarkerSize',9, ...
             'DisplayName','periodic orbit (fixed point)')
    else
        plot((r - r0)*c.DU, rd*c.VU*1000, '.', 'MarkerSize',6, ...
             'DisplayName',sprintf('start + %g m/s', kick(k)*c.VU*1000))
    end
    fprintf('  %s Poincare map, kick %4.1f m/s: r spread %.3e km, rdot spread %.3e m/s\n', ...
            name, kick(k)*c.VU*1000, (max(r)-min(r))*c.DU, (max(rd)-min(rd))*c.VU*1000);
end
xlabel('r - r_0 at ascending node  [km]'); ylabel('radial speed at node  [m/s]')
title(sprintf('%s: Poincare map, %d laps', name, nlap))
legend('Location','eastoutside')
fig_format(gcf, 24)

end


function [val,stop,dir] = event_node_all(~,X)
val  = X(3);    
stop = 0;      
dir  = 1;      
end

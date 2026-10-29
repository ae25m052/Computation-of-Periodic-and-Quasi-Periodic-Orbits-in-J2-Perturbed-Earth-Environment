% Time from one perigee to the next for 100 revolutions,
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
T  = 2*pi*sqrt(a^3/mu);           % Kepler period, TU

nrev = 100;

i  = acos(sqrt(1/5));             % critical inclination from the equation
Wf = -1.5*n*J2*(Re/p)^2*cos(i);   % frame rate = nodal rate

% start at perigee on the +x axis
[r0,v0] = kep2rv(a,e,i,0,0,0,mu);

opts = odeset(c.opts,'Events',@event_peri_dir);
tend = (nrev+1)*T;                % a bit more than 100 revolutions

% ECI
[~,~,te,Xe,ie] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,1),[0 tend],[r0;v0],opts);
[tpI,tdI] = split_events(te,Xe,ie,nrev);

% rotating frame, same start, velocity seen from the frame
X0r = [r0; v0 - cross([0;0;Wf],r0)];
[~,~,te,Xe,ie] = ode89(@(t,X) rf_eom(t,X,mu,Re,J2,Wf),[0 tend],X0r,opts);
[tpR,tdR] = split_events(te,Xe,ie,nrev);

% per revolution: back to start direction, perigee to perigee (min),
% extra time in that revolution and total extra up to it (s)
dI = diff(tdI)*c.TU/60;   PI = diff(tpI)*c.TU/60;
dR = diff(tdR)*c.TU/60;   PR = diff(tpR)*c.TU/60;
xI = (PI - dI)*60;        sI = cumsum(xI);
xR = (PR - dR)*60;        sR = cumsum(xR);

fprintf('hp = %.1f km, e = %.2f, i = %.6f deg, Wf = %+.4f deg/day\n', ...
        hp,e,rad2deg(i),rad2deg(Wf/c.TU*86400));
fprintf('Kepler period T = %.6f min\n',T*c.TU/60);
fprintf('back to start = time to come back to the starting direction (y = 0, x > 0)\n');
fprintf('extra         = perigee to perigee - back to start (extra time in that revolution)\n');
fprintf('total extra   = sum of extra up to that revolution\n');
fprintf('rev 0         = starting point, at perigee, t = 0\n\n');
fprintf('         ------------------------ ECI ------------------------     ------------------- rotating frame --------------------\n');
fprintf('  rev   back to start  extra    total extra  perigee to perigee    back to start  extra    total extra  perigee to perigee\n');
fprintf('          [min]         [s]        [s]          [min]                [min]         [s]        [s]          [min]\n');
fprintf('  %3d   %10s   %7s   %+9.4f   %12s       %10s   %7s   %+9.4f   %12s\n', ...
        0,'-','-',0,'-','-','-',0,'-');
for k = 1:nrev
    fprintf('  %3d   %.6f   %+.4f   %+9.4f   %.8f       %.6f   %+.4f   %+9.4f   %.8f\n', ...
            k,dI(k),xI(k),sI(k),PI(k),dR(k),xR(k),sR(k),PR(k));
end

fprintf('\nperigee to perigee, max - min over %d revolutions: ECI %.3e s, rotating %.3e s\n', ...
        nrev,(max(PI)-min(PI))*60,(max(PR)-min(PR))*60);
fprintf('max |ECI - rotating| perigee to perigee = %.3e s\n',max(abs(PI-PR))*60);
fprintf('mean perigee to perigee - Kepler T      = %+.4f s\n',(mean(PI) - T*c.TU/60)*60);

figure
subplot(2,1,1)
plot(1:nrev,PI,'b'); hold on; grid on
plot(1:nrev,PR,'r--');
ylabel('perigee to perigee [min]')
legend('ECI','rotating frame','Location','best')
subplot(2,1,2)
plot(1:nrev,xI,'b'); hold on; grid on
plot(1:nrev,xR,'r--');
xlabel('revolution'); ylabel('extra time per revolution [s]')
legend('ECI','rotating frame','Location','best')
fig_format(gcf)


function [value,isterminal,direction] = event_peri_dir(~,X)
% 1: perigee, r.v goes from negative to positive
% 2: starting direction, y goes from negative to positive (x > 0 checked later)

value      = [X(1:3).'*X(4:6); X(2)];
isterminal = [0; 0];
direction  = [1; 1];

end


function [tp,td] = split_events(te,Xe,ie,nrev)
% Perigee times and start-direction times, with t = 0 as number 0.

k  = te > 1e-6;
te = te(k);  Xe = Xe(k,:);  ie = ie(k);

tp = [0; te(ie == 1)];
td = [0; te(ie == 2 & Xe(:,1) > 0)];

tp = tp(1:nrev+1);
td = td(1:nrev+1);

end

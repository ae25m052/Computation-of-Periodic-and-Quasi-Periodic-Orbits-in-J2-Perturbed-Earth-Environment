% Rotating frame check. The same orbit is solved in ECI (j2_eom) and in a
% frame spinning about z at rate Wf (rf_eom). 
clear; clc; close all

c = j2_constants();

% normalise first
mu = 1;
Re = 1;
J2 = c.J2;

hp = 300;                         % km
e  = 0.65;
i  = deg2rad(63.43);
rp = (c.Re + hp)/c.DU;            % DU
a  = rp/(1-e);
p  = a*(1-e^2);
n  = sqrt(mu/a^3);
T  = 2*pi*sqrt(a^3/mu);           % TU

Wf = -1.5*n*J2*(Re/p)^2*cos(i);   % frame rate = nodal rate

nrev = 50;
t = linspace(0,nrev*T,2001);

[r0,v0] = kep2rv(a,e,i,0,0,0,mu);

% ECI
[~,XI] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,1),t,[r0;v0],c.opts);

% rotating frame: same position, velocity seen from the frame
W   = [0; 0; Wf];
X0r = [r0; v0 - cross(W,r0)];
[~,XR] = ode89(@(t,X) rf_eom(t,X,mu,Re,J2,Wf),t,X0r,c.opts);

% rotate back to ECI:  r = Rz(Wf t) * rho
rb = zeros(numel(t),3);
for k = 1:numel(t)
    th = Wf*t(k);
    R  = [cos(th) -sin(th) 0
          sin(th)  cos(th) 0
          0        0       1];
    rb(k,:) = (R*XR(k,1:3).').';
end
err = sqrt(sum((rb - XI(:,1:3)).^2,2));

% Jacobi constant in the rotating frame
x = XR(:,1);  y = XR(:,2);  z = XR(:,3);
r = sqrt(x.^2 + y.^2 + z.^2);
V = -mu./r + mu*J2*Re^2./(2*r.^3).*(3*z.^2./r.^2 - 1);
C = 0.5*sum(XR(:,4:6).^2,2) + V - 0.5*Wf^2*(x.^2 + y.^2);

fprintf('hp = %.1f km, e = %.2f, i = %.2f deg, %d revolutions\n',hp,e,rad2deg(i),nrev);
fprintf('Wf = %+.6e rad/TU = %+.4f deg/day\n\n',Wf,rad2deg(Wf/c.TU*86400));
fprintf('max |r_ECI - Rz*rho| = %.3e DU = %.3e km\n',max(err),max(err)*c.DU);
fprintf('Jacobi constant C(0) = %.12f\n',C(1));
fprintf('max |C - C(0)|       = %.3e\n',max(abs(C - C(1))));

figure
subplot(2,1,1)
semilogy(t/T,err*c.DU,'b'); grid on
ylabel('position difference [km]')
subplot(2,1,2)
plot(t/T,C - C(1),'b'); grid on
xlabel('revolution'); ylabel('C - C(0)')
fig_format(gcf)

plot_orbit(t,XR,T,'Rotating frame');
fig_format(gcf)
plot_orbit(t,XI,T,'ECI');
fig_format(gcf)

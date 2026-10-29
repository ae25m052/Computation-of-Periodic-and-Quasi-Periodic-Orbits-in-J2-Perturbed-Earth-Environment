function [vp,i,Wf,T,hist] = shoot_inclined(rp,ra,i,mu,Re,J2,opts)
% Fixed i (ra = []): unknowns vp, Wf      conditions xdot = 0, zdot = 0
% Free i  (i = [])  : unknowns vp, i, Wf   conditions xdot = 0, zdot = 0, |r| - ra = 0

free_i = isempty(i);

if free_i
    e = (ra-rp)/(ra+rp);
    i = deg2rad(63.43);
else
    e = 0;
end

% Kepler speed and Brouwer nodal rate as first guess
a  = rp/(1-e);
p  = a*(1-e^2);
n  = sqrt(mu/a^3);
vp = sqrt(mu*(1+e)/rp);
Wf = -1.5*n*J2*(Re/p)^2*cos(i);

if free_i
    u = [vp; i; Wf];
else
    u = [vp; Wf];
end
m = numel(u);

hist = [];

for k = 1:20
    [F,th] = half_lap(u,rp,ra,i,free_i,mu,Re,J2,opts);
    hist = [hist; k u.' max(abs(F))];

    if max(abs(F)) < 1e-10
        break
    end

    J = zeros(m);
    for j = 1:m
        h  = 1e-7*max(abs(u(j)),1e-3);
        uh = u;  uh(j) = uh(j) + h;
        J(:,j) = (half_lap(uh,rp,ra,i,free_i,mu,Re,J2,opts) - F)/h;
    end

    u = u - J\F;
end

vp = u(1);
Wf = u(end);
if free_i
    i = u(2);
end
T = 2*th;

end


function [F,th] = half_lap(u,rp,ra,i,free_i,mu,Re,J2,opts)
% Propagate in the rotating frame until y = 0 again and return the
% conditions F there and the time th (half the period).
vp = u(1);
Wf = u(end);
if free_i
    i = u(2);
end

X0 = [rp*cos(i); 0; rp*sin(i); 0; vp - Wf*rp*cos(i); 0];

a_guess = 1/(2/rp - vp^2/mu);
tmax    = 2*pi*sqrt(a_guess^3/mu);

opts = odeset(opts,'Events',@event_ycross);
[~,~,te,Xe] = ode89(@(t,X) rf_eom(t,X,mu,Re,J2,Wf),[0 tmax],X0,opts);

th = te(1);
Xh = Xe(1,:);

F = [Xh(4); Xh(6)];
if free_i
    F = [F; norm(Xh(1:3)) - ra];
end

end


function [value,isterminal,direction] = event_ycross(~,X)
% Orbit back at the x-z plane (y = 0, y decreasing).
value      = X(2);
isterminal = 1;
direction  = -1;

end

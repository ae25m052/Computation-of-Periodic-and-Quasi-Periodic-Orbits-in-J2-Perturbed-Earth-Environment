function [ra,thalf] = apogee_radius(vp,rp,mu,Re,J2,opts)
% Start at perigee on the +x axis with speed vp and propagate until r.v
% changes sign, which is apoapsis. Returns the apogee radius and the time
% taken to get there (half the radial period).

X0 = [rp; 0; 0; 0; vp; 0];

a_guess = 1/(2/rp - vp^2/mu);
tmax    = 2*pi*sqrt(a_guess^3/mu);

opts = odeset(opts,'Events',@event_apsis);
[~,~,te,Xe] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,1),[0 tmax],X0,opts);

if isempty(te)
    ra    = NaN;
    thalf = NaN;
else
    ra    = norm(Xe(1,1:3));
    thalf = te(1);
end

end

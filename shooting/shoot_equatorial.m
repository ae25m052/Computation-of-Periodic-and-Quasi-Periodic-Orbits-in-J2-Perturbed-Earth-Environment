function [vp,ra,Trad,hist] = shoot_equatorial(rp,ra_target,mu,Re,J2,opts)
% Newton correction on the perigee speed so that the J2 orbit reaches the
% intended apogee radius. Perigee radius is held fixed.

e_guess = (ra_target-rp)/(ra_target+rp);
vp      = sqrt(mu*(1+e_guess)/rp);          % Kepler starting guess

hist = [];

for k = 1:20
    [ra,thalf] = apogee_radius(vp,rp,mu,Re,J2,opts);
    fval = ra - ra_target;
    hist = [hist; k vp fval];

    if abs(fval) < 1e-10
        break
    end

    h    = 1e-7*vp;
    ra_h = apogee_radius(vp+h,rp,mu,Re,J2,opts);
    dfdv = (ra_h - ra)/h;

    vp = vp - fval/dfdv;
end

Trad = 2*thalf;

end

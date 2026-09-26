%MAIN_STAGE2B_ANALYTIC_STM  Can the STM be derived analytically?
%   THREE DIFFERENT THINGS "ANALYTICAL STM" CAN MEAN
%
%     (a) An analytical JACOBIAN A(t).
%         Already done in Stage 2 -- gravity_gradient.m is a hand-derived
%         closed form, exactly symmetric.  Nothing is differenced.
%
%     (b) Solving Phi' = A(t) Phi symbolically, i.e. Phi = exp(int A dt).
%         FAILS.  Section 2 below measures the failure: 665 % error even for
%         pure Kepler motion, because A(t) does not commute with itself at
%         different times.
%
%     (c) Differentiating a known closed-form solution of the NONLINEAR flow.
%         WORKS whenever that flow is integrable in closed form.  Section 1
%         does exactly this for Kepler and matches numerical propagation to
%         ~1e-13 relative.  Section 3 shows why it has no J2 counterpart.

clear; clc;

mu_dim = 398600.4418;  Re_dim = 6378.1363;  J2 = 1.08262668e-3;
DU = Re_dim;  TU = sqrt(DU^3/mu_dim);
mu = 1;  Re = 1;
opts = odeset('RelTol',1e-13,'AbsTol',1e-13);

oe0 = [7500/DU; 0.05; deg2rad(45); deg2rad(30); deg2rad(20); 0];
T   = 2*pi*sqrt(oe0(1)^3/mu);

% initial Cartesian state from the same element set
J0 = dxdoe(oe0, mu);   %#ok<NASGU>  (also used below)
x0 = oe2rv_local(oe0, mu);

fprintf('=== STAGE 2b: ANALYTICAL STM INVESTIGATION ===\n\n');

%% ---------------------------------------------------------------
%  0. The analytical element partials are correct
%  -----------------------------------------------------------------
Ja = dxdoe(oe0, mu);
Jf = zeros(6,6);  h = 1e-7;
for k = 1:6
    op = oe0; om = oe0;  op(k) = op(k)+h;  om(k) = om(k)-h;
    Jf(:,k) = (oe2rv_local(op,mu) - oe2rv_local(om,mu))/(2*h);
end
fprintf('--- 0. Analytical element partials d(r,v)/d(oe) ---\n');
fprintf('  max |analytical - finite difference| = %.3e   (FD floor at h=1e-7)\n\n', ...
        max(abs(Ja-Jf),[],'all'));

%% ---------------------------------------------------------------
%  1. ANALYTICAL KEPLERIAN STM -- this works, exactly
%  -----------------------------------------------------------------
fprintf('--- 1. Closed-form Keplerian STM vs numerical propagation ---\n');
for k = [1 5 20]
    [Pa, N] = stm_kepler_analytic(oe0, k*T, mu);
    Pn = stm_numeric(x0, k*T, 0, mu, J2, Re, opts);
    fprintf('  %2d orbits : max|dPhi| = %.3e ,  relative = %.3e\n', ...
            k, max(abs(Pa-Pn),[],'all'), ...
            max(abs(Pa-Pn),[],'all')/max(abs(Pn),[],'all'));
end

[Pa, N] = stm_kepler_analytic(oe0, T, mu);
fprintf('\n  Structure, derived analytically rather than observed:\n');
fprintf('    N is nilpotent, max|N*N|     = %.1e   (exactly zero)\n', max(abs(N*N),[],'all'));
fprintf('    rank(Phi - I)                = %d\n', rank(Pa-eye(6),1e-8));
fprintf('    max|(Phi - I)^2|             = %.3e\n', max(abs((Pa-eye(6))^2),[],'all'));
fprintf('    -> Jordan block at +1. This EXPLAINS the 821x eigenvalue split\n');
fprintf('       measured in Stage 2: defective eigenvalues split like sqrt(eps).\n\n');

%% ---------------------------------------------------------------
%  2. THE OBVIOUS SHORTCUT, AND WHY IT FAILS
%  -----------------------------------------------------------------
% Phi' = A(t) Phi looks like a matrix version of y' = a(t) y, whose solution
% is exp(int a dt).  That step is legal ONLY if A(t1) and A(t2) commute for
% every pair of times.  They do not, and the damage is not subtle.
fprintf('--- 2. Does Phi = exp( int A dt ) work? ---\n');
for eps = [0 1]
    [tt, xx] = ode89(@(t,x) twobody_eom(t,x,mu,eps,J2,Re), linspace(0,T,4001), x0, opts);

    ng = 9;  ii = round(linspace(1,numel(tt),ng));
    cmax = 0;  amax = 0;
    Ai = cell(1,ng);
    for p = 1:ng
        Ai{p} = Amat(xx(ii(p),1:3).', mu, eps, J2, Re);
        amax  = max(amax, max(abs(Ai{p}),[],'all'));
    end
    for p = 1:ng
        for q = 1:ng
            cmax = max(cmax, max(abs(Ai{p}*Ai{q} - Ai{q}*Ai{p}),[],'all'));
        end
    end

    Aint = zeros(6,6);
    Astack = zeros(numel(tt),36);
    for k = 1:numel(tt)
        Ak = Amat(xx(k,1:3).', mu, eps, J2, Re);
        Astack(k,:) = Ak(:).';
    end
    Aint(:) = trapz(tt, Astack, 1);

    Pe = expm(Aint);
    Pn = stm_numeric(x0, T, eps, mu, J2, Re, opts);

    fprintf('  eps = %d :\n', eps);
    fprintf('    max |[A(t1), A(t2)]|            = %.4e   (would be 0 if exp were valid)\n', cmax);
    fprintf('    max |A|                         = %.4e\n', amax);
    fprintf('    commutator / |A|^2              = %.4f   (order one, not small)\n', cmax/amax^2);
    fprintf('    max |exp(int A dt) - Phi_true|  = %.4e\n', max(abs(Pe-Pn),[],'all'));
    fprintf('    RELATIVE ERROR                  = %.1f %%\n', ...
            100*max(abs(Pe-Pn),[],'all')/max(abs(Pn),[],'all'));
end
fprintf('\n  Note this fails at eps = 0 too, where a closed-form STM DOES exist.\n');
fprintf('  So the obstruction is not J2: it is that exponentiating a time-varying\n');
fprintf('  matrix requires the Magnus series, which does not terminate here.\n\n');

%% ---------------------------------------------------------------
%  3. WHY THERE IS NO J2 COUNTERPART TO SECTION 1
%  -----------------------------------------------------------------
% Section 1 worked because the Keplerian flow is closed-form in element space:
% five elements are constant and the sixth advances linearly.  Under J2 the
% osculating elements do neither -- they oscillate within every revolution.
fprintf('--- 3. Osculating elements under J2, over ONE revolution ---\n');
[tt, xx] = ode89(@(t,x) twobody_eom(t,x,mu,1,J2,Re), linspace(0,T,400), x0, opts);
av = zeros(size(tt)); ev = av; iv = av;
for k = 1:numel(tt)
    r = xx(k,1:3).';  v = xx(k,4:6).';  rn = norm(r);
    hv = cross(r,v);
    av(k) = 1/(2/rn - dot(v,v)/mu);
    ev(k) = norm(cross(v,hv)/mu - r/rn);
    iv(k) = acos(hv(3)/norm(hv));
end
fprintf('  a varies by %.3f km,  e by %.3e,  i by %.5f deg\n', ...
        (max(av)-min(av))*DU, max(ev)-min(ev), rad2deg(max(iv)-min(iv)));
fprintf('  These are the short-period J2 terms. There is no closed-form\n');
fprintf('  expression for them: the main problem of artificial satellite theory\n');
fprintf('  is not integrable. Brouwer''s theory gives MEAN elements in closed\n');
fprintf('  form to first order in J2, which is what an approximate analytical\n');
fprintf('  STM (Gim & Alfriend 2003) is built on -- and its error was measured\n');
fprintf('  in Stage 2 as 1.2%% after one orbit rising to 18.3%% after twenty.\n');

%% ---------------------------------------------------------------
%  Local functions
%  -----------------------------------------------------------------
function Phi = stm_numeric(x0, tf, eps, mu, J2, Re, opts)
Y0 = [x0; reshape(eye(6),36,1)];
[~, Y] = ode89(@(t,Y) stm_eom(t,Y,mu,eps,J2,Re), [0 tf], Y0, opts);
Phi = reshape(Y(end,7:42), 6, 6);
end

function A = Amat(r, mu, eps, J2, Re)
G = gravity_gradient(r, mu, eps, J2, Re);
A = [zeros(3) eye(3) ; G zeros(3)];
end

function x = oe2rv_local(oe, mu)
a=oe(1); e=oe(2); i=oe(3); Om=oe(4); w=oe(5); M=oe(6);
E = M;  if e > 0.8, E = pi; end
for k = 1:100
    d = (E - e*sin(E) - M)/(1 - e*cos(E));  E = E - d;
    if abs(d) < 1e-15, break; end
end
s=sin(E); c=cos(E); b=sqrt(1-e^2);
r=a*(1-e*c); eta=sqrt(mu*a);
p=[a*(c-e); a*b*s; 0];
q=[-eta*s/r; eta*b*c/r; 0];
cO=cos(Om); sO=sin(Om); ci=cos(i); si=sin(i); cw=cos(w); sw=sin(w);
R=[ cO*cw-sO*sw*ci , -cO*sw-sO*cw*ci ,  sO*si ; ...
    sO*cw+cO*sw*ci , -sO*sw+cO*cw*ci , -cO*si ; ...
    sw*si          ,  cw*si          ,  ci    ];
x = [R*p; R*q];
end

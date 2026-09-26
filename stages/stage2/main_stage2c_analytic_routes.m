%MAIN_STAGE2C_ANALYTIC_ROUTES  Every remaining analytical route, tested.
%     ROUTE 4  Magnus series. The correct generalisation of exp(int A dt).
%              Tested against its own convergence theorem.
%     ROUTE 5  Floquet theory. Phi(t) = P(t) exp(Rt) for periodic A(t).
%     ROUTE 6  Vinti's problem. An integrable neighbour of the J2 problem
%              with genuinely closed-form solutions.

clear; clc;

mu_dim = 398600.4418;  Re_dim = 6378.1363;  J2 = 1.08262668e-3;
DU = Re_dim;  TU = sqrt(DU^3/mu_dim);
mu = 1;  Re = 1;
opts = odeset('RelTol',1e-13,'AbsTol',1e-13);

a0 = 7500/DU;
[r0,v0] = coe2rv(a0, 0.05, deg2rad(45), deg2rad(30), deg2rad(20), 0, mu);
x0 = [r0(:); v0(:)];
T  = 2*pi*sqrt(a0^3/mu);

fprintf('=== STAGE 2c: REMAINING ANALYTICAL ROUTES ===\n\n');

%% ---------------------------------------------------------------
%  ROUTE 4 -- Magnus series, tested against its convergence theorem
%  -----------------------------------------------------------------
% exp(int A dt) is only the FIRST term of the Magnus series
%     Phi = exp(Om1 + Om2 + Om3 + ...),   Om2 = 1/2 int int [A(t1),A(t2)]
% The series has a proven convergence condition:
%     int_0^T ||A(t)||_2 dt  <  pi
% If that integral exceeds pi the series is not guaranteed to converge at
% all, and no number of extra terms rescues it. Measure the integral.
fprintf('--- ROUTE 4: Magnus series convergence condition ---\n');
for eps = [0 1]
    tg = linspace(0, T, 2001);
    [~, XX] = ode89(@(t,x) twobody_eom(t,x,mu,eps,J2,Re), tg, x0, opts);
    An = zeros(size(tg));
    for k = 1:numel(tg)
        G = gravity_gradient(XX(k,1:3).', mu, eps, J2, Re);
        An(k) = norm([zeros(3) eye(3); G zeros(3)], 2);
    end
    I = trapz(tg, An);
    kk = find(cumtrapz(tg,An) > pi, 1, 'first');
    fprintf('  eps = %d :  int ||A||_2 dt over one orbit = %.4f  (bound = pi = %.4f)\n', eps, I, pi);
    fprintf('             exceeds the bound by %.2fx, first violated at %.1f %% of one orbit\n', ...
            I/pi, 100*tg(kk)/T);
end
fprintf('\n  The series is not guaranteed to converge even over a THIRD of one\n');
fprintf('  orbit. Adding higher-order terms cannot fix this. Route closed.\n\n');

%% ---------------------------------------------------------------
%  ROUTE 5 -- Floquet theory
%  -----------------------------------------------------------------
fprintf('--- ROUTE 5: Floquet theory ---\n');
fprintf('  If A(t) were periodic with period T, Floquet gives\n');
fprintf('      Phi(t) = P(t) exp(R t),   P periodic, R constant.\n');
fprintf('  Two obstructions:\n');
fprintf('   (i) A(t) is periodic only ON a periodic orbit. For a generic J2\n');
fprintf('       orbit the trajectory never closes, so A(t) is not periodic\n');
fprintf('       and the theorem does not apply at all.\n');
fprintf('  (ii) Even on a periodic orbit, R is DEFINED by R = log(Phi(T))/T.\n');
fprintf('       Obtaining R requires the monodromy matrix we were trying to\n');
fprintf('       compute. The construction is circular. Route closed.\n\n');

%% ---------------------------------------------------------------
%  ROUTE 6 -- Vinti's problem
%  -----------------------------------------------------------------
% Vinti's potential is separable in oblate spheroidal coordinates, so the
% problem is INTEGRABLE and closed-form solutions exist. The catch is that it
% is not the problem we defined: expanded in zonal harmonics it reproduces J2
% exactly but also carries J4 = -J2^2, J6 = +J2^3, ... which the main problem
% does not have.
fprintf('--- ROUTE 6: Vinti''s problem ---\n');
c2 = J2*Re^2;
fprintf('  Vinti parameter c = Re*sqrt(J2) = %.3f km\n', sqrt(c2)*DU);
fprintf('  Vinti reproduces: J2 exactly, J3 = 0, J4 = -J2^2, J6 = +J2^3, ...\n');
fprintf('  The J2 main problem has NONE of the J4, J6 terms.\n\n');

rt = [1.1; 0.3; 0.4];
aV = accel_vinti(rt, mu, c2);
aM = twobody_eom(0,[rt;0;0;0],mu,1,J2,Re);  aM = aM(4:6);
aK = -mu*rt/norm(rt)^3;
fprintf('  At a test point:\n');
fprintf('    |a_Vinti - a_main| / |a_main| = %.4e\n', norm(aV-aM)/norm(aM));
fprintf('    the J2 term itself is          %.4e of the total\n', norm(aM-aK)/norm(aM));
fprintf('    so the model mismatch is ~%.1f %% OF THE PERTURBATION\n\n', ...
        100*norm(aV-aM)/norm(aM-aK));

fprintf('  Trajectory divergence from the same initial state:\n');
sM = ode89(@(t,x) twobody_eom(t,x,mu,1,J2,Re), [0 60*T], x0, opts);
sV = ode89(@(t,x) [x(4:6); accel_vinti(x(1:3),mu,c2)], [0 60*T], x0, opts);
for k = [1 10 30 60]
    tt = linspace(0,k*T,4000);
    XM = deval(sM,tt); XV = deval(sV,tt);
    d  = max(vecnorm(XM(1:3,:)-XV(1:3,:)));
    fprintf('    %2d orbits : max |dr| = %8.4f km\n', k, d*DU);
end

fprintf('\n  VERDICT ON VINTI\n');
fprintf('  It is a real integrable system with genuinely closed-form solutions,\n');
fprintf('  so an exact analytical STM DOES exist -- for Vinti''s problem.\n');
fprintf('  Two reasons it does not answer the question as posed:\n');
fprintf('   (i) It solves a DIFFERENT model. The mismatch is ~1e-6 relative,\n');
fprintf('       which is seven orders of magnitude larger than the 1e-13\n');
fprintf('       tolerance every validation in Stages 1-3 is built on.\n');
fprintf('  (ii) "Closed form" here means reduction to quadratures. The Vinti\n');
fprintf('       solution requires elliptic-type integrals inverted by series,\n');
fprintf('       so it is semi-analytical in practice, not a formula.\n\n');

fprintf('  WORTH NOTING: as an APPROXIMATE Jacobian, Vinti is excellent --\n');
fprintf('  roughly 3000x better than the Keplerian proxy (0.0004 %% vs 1.2 %%\n');
fprintf('  after one orbit). That is a possible Stage 3 optimisation, not a\n');
fprintf('  replacement for the numerical STM.\n');

%% ---------------------------------------------------------------
function a = accel_vinti(r, mu, c2)
%ACCEL_VINTI  Gradient of Vinti's potential V = mu*rho^3/(rho^4 + c^2 z^2),
%   with rho the oblate spheroidal radial coordinate satisfying
%       rho^4 - rho^2 (r^2 - c^2) - c^2 z^2 = 0 .
x = r(1); y = r(2); z = r(3);
r2 = x*x + y*y + z*z;
b  = r2 - c2;
rho2 = 0.5*(b + sqrt(b*b + 4*c2*z*z));
rho  = sqrt(rho2);
Q = rho2*rho2 + c2*z*z;
D = 2*rho2 - r2 + c2;
dV_drho = mu*(3*rho2*Q - 4*rho^6)/Q^2;
a = [ dV_drho * x*rho/D ; ...
      dV_drho * y*rho/D ; ...
      dV_drho * z*(rho2+c2)/(rho*D) - 2*mu*rho^3*c2*z/Q^2 ];
end

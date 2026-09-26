function [Phi, N] = stm_kepler_analytic(oe0, t, mu)
%STM_KEPLER_ANALYTIC  Closed-form Keplerian STM. No integration whatsoever.
%
%   [Phi, N] = stm_kepler_analytic(oe0, t, mu)
%
%   oe0 = [a; e; i; RAAN; argp; M0],  t = elapsed time.
%
%   DERIVATION
%   In Cartesian coordinates the Keplerian flow is not obviously integrable,
%   but in classical elements it is trivial:
%
%       a, e, i, RAAN, argp   all constant
%       M(t) = M0 + n(a) t ,   n(a) = sqrt(mu/a^3)
%
%   so the STM in ELEMENT space is exactly
%
%       Phi_oe(t) = I + N t ,    N(M,a) = dn/da = -3n/(2a) ,
%
%   a single off-diagonal entry.  Transforming back through the analytical
%   element partials of dxdoe.m,
%
%       Phi(t) = J(oe(t)) * (I + N t) * J(oe(0))^{-1} .
%
%   WHY THE MONODROMY IS DEFECTIVE -- the analytical explanation
%   N is NILPOTENT: N^2 = 0 exactly, because its only nonzero entry is
%   off-diagonal.  So Phi_oe = I + Nt is a Jordan block, not a diagonalizable
%   matrix, and all six eigenvalues sit at +1 with only ONE independent
%   eigenvector.  Similarity by J preserves this, giving
%
%       rank(Phi - I) = 1 ,   (Phi - I)^2 = 0 ,
%
%   exactly as measured numerically in Stage 2.  A defective eigenvalue splits
%   under perturbation like sqrt(eps) rather than linearly, which is why a
%   perturbation of order 1e-3 produces the observed factor-821 separation.
%   That is not a numerical artefact; it is the structure of the problem.
%
%   VALIDITY
%   Exact for eps = 0 only.  There is no counterpart for eps = 1: see
%   main_stage2b_analytic_stm.m.

a = oe0(1);
n = sqrt(mu/a^3);

oet    = oe0;
oet(6) = oe0(6) + n*t;              % only the mean anomaly advances

N = zeros(6,6);
N(6,1) = -1.5*n/a;                  % dn/da

Phi = dxdoe(oet, mu) * (eye(6) + N*t) / dxdoe(oe0, mu);
end

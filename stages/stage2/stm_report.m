function R = stm_report(Phi)
%STM_REPORT  Everything worth measuring about one state transition matrix.
%
%   R = stm_report(Phi)
%
%   Output struct:
%       ev            6x1 eigenvalues, sorted by magnitude (ascending)
%       absev         6x1 magnitudes, ascending
%       maxabs        largest magnitude
%       symp          max |Phi'*J*Phi - J|
%       det1          |det(Phi) - 1|
%       pair_err      worst reciprocal-pair error over the three pairs
%       n_on_circle   how many eigenvalues lie on the unit circle to 1e-6
%       rank_PhiI     rank(Phi - I)
%       jordan_resid  ||(Phi-I)^2|| / ||Phi-I||^2
%
%   ON SYMPLECTICITY -- the sharpest test available
%   A Hamiltonian flow satisfies Phi'*J*Phi = J exactly, so whatever this
%   prints is pure integrator error and nothing else. It is far sharper
%   than comparing against finite differences (which bottoms out near 1e-8
%   because of its own step size), and it is a JOINT test of the
%   acceleration and its Jacobian: a single wrong coefficient in the
%   analytically derived A matrix leaves the trajectory looking perfectly
%   plausible while this number degrades visibly.
%
%   Expect ~1e-11 for the e = 0.65, h_p = 300 km orbit. That is NOT a
%   regression from the ~1e-13 seen on the earlier near-circular 7500 km
%   orbit -- it is the cost of a fast perigee passage at high eccentricity,
%   and it does not improve below RelTol = 1e-13.
%
%   ON RECIPROCAL PAIRING
%   Symplectic matrices have eigenvalues in pairs (lambda, 1/lambda).
%   Sorting by magnitude and multiplying the smallest by the largest, the
%   second by the fifth, and so on must give 1 three times over. Unlike the
%   "all eigenvalues equal 1" test, this one keeps working at epsilon = 1,
%   where the eigenvalues have separated. It is the structural check that
%   survives into the perturbed problem.
%
%   ON rank_PhiI AND jordan_resid -- THE RIGHT WAY TO TEST A KEPLERIAN STM
%   At epsilon = 0 the STM is a SHEAR: a single Jordan block, all six
%   eigenvalues exactly 1, but only one independent eigenvector. The two
%   structural signatures of that are
%
%       rank(Phi - I) = 1            and         (Phi - I)^2 = 0
%
%   Both are computed here, the second normalised by ||Phi-I||^2 so it
%   reads as a relative error.
%
%   These are more trustworthy than checking that the eigenvalues equal 1,
%   because computing the eigenvalues of a defective matrix is an
%   ill-conditioned problem BY NATURE: the eigenvalues of a Jordan block
%   respond to a perturbation of size d like d^(1/2), not like d. So a
%   Keplerian STM whose entries are correct to 1e-13 can still report
%   eigenvalues that miss 1 by something in the 1e-6 range, and different
%   eigenvalue routines on the same matrix can disagree at that level.
%   Nothing is wrong when that happens -- the eigenvalue solver is simply
%   being asked a question the matrix cannot answer sharply. The rank and
%   nilpotency tests go nowhere near that difficulty, so they are what
%   should be believed.
%
%   M.Tech. project, Raj Khismatrao (AE25M052).

    Jsym = [zeros(3,3), eye(3);
            -eye(3),    zeros(3,3)];

    ev = eig(Phi);
    [absev, idx] = sort(abs(ev));

    R.ev          = ev(idx);
    R.absev       = absev;
    R.maxabs      = absev(end);
    R.symp        = max(max(abs(Phi.' * Jsym * Phi - Jsym)));
    R.det1        = abs(det(Phi) - 1);
    R.n_on_circle = sum(abs(absev - 1) < 1e-6);

    pe = zeros(3,1);
    for m = 1:3
        pe(m) = abs(abs(R.ev(m) * R.ev(7-m)) - 1);
    end
    R.pair_err = max(pe);

    % ---- structural (Jordan) diagnostics --------------------------------
    M  = Phi - eye(6);
    nM = norm(M);
    if nM > 0
        R.rank_PhiI    = rank(M, 1e-6*nM);
        R.jordan_resid = norm(M*M) / nM^2;
    else
        R.rank_PhiI    = 0;
        R.jordan_resid = 0;
    end
end

function [energy, h_vec, e_vec] = compute_invariants(x, mu)
%   x   - state vector [r(3); v(3)], one column per time sample (6xN)
%   mu  - gravitational parameter (matching units/non-dim of x)
%
%   Outputs (each one column per time sample):
%       energy - specific orbital energy,  v^2/2 - mu/r        (1xN)
%       h_vec  - specific angular momentum vector, r x v        (3xN)
%       e_vec  - eccentricity (LRL) vector, (v x h)/mu - r/|r|  (3xN)

    N = size(x, 2);
    energy = zeros(1, N);
    h_vec  = zeros(3, N);
    e_vec  = zeros(3, N);

    for k = 1:N
        r_vec = x(1:3, k);
        v_vec = x(4:6, k);
        r = norm(r_vec);
        v = norm(v_vec);

        energy(k) = v^2/2 - mu/r;

        h = cross(r_vec, v_vec);
        h_vec(:,k) = h;

        e_vec(:,k) = cross(v_vec, h)/mu - r_vec/r;
    end
end

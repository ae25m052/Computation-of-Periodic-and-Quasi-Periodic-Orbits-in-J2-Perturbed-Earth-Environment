function H = element_history(X, mu)
%ELEMENT_HISTORY  Osculating elements along a state history, angles unwrapped.
%
%   H = element_history(X, mu)      X is 6xN
%
%   Output struct of 1xN arrays plus scalars:
%       a            semi-major axis, same units as X
%       e            eccentricity
%       inc_deg      inclination                          [deg]
%       RAAN_deg     right ascension of node, wrapped     [deg]
%       argp_deg     argument of perigee, wrapped         [deg]
%       RAAN_un_deg  RAAN, UNWRAPPED                      [deg]
%       argp_un_deg  argp, UNWRAPPED                      [deg]
%       d_RAAN_deg   total drift, last minus first        [deg]
%       d_argp_deg   total drift, last minus first        [deg]
%       equatorial   TRUE if the equatorial branch was used anywhere
%
%   WHY UNWRAPPING MATTERS
%   A wrapped angle jumps from 360 back to 0 every time it goes round, so a
%   raw plot of it is a sawtooth and the total drift cannot be read off at
%   all. unwrap() removes those jumps, turning the sawtooth into the
%   straight line whose slope is the secular rate. That line is the whole
%   point of the exercise: it is how the critical inclination announces
%   itself, as the one inclination where the argp line is flat.
%
%   WHAT TO EXPECT UNDER J2
%   a, e and i have NO secular drift -- they only oscillate with short
%   period. RAAN and argp do drift secularly. That contrast is itself a
%   check on the propagation.
%
%   If H.equatorial is TRUE the orbit was equatorial somewhere along the
%   history, in which case argp is really the longitude of perigee. See
%   rv2coe_safe.m for why, and label plots accordingly.
%
%   M.Tech. project, Raj Khismatrao (AE25M052).

    N = size(X, 2);

    a_ = zeros(1,N);  e_ = zeros(1,N);  i_ = zeros(1,N);
    O_ = zeros(1,N);  w_ = zeros(1,N);
    eqf = false;

    for k = 1:N
        [a_(k), e_(k), i_(k), O_(k), w_(k), ~, eq] = ...
            rv2coe_safe(X(1:3,k), X(4:6,k), mu);
        eqf = eqf || eq;
    end

    H.a           = a_;
    H.e           = e_;
    H.inc_deg     = rad2deg(i_);
    H.RAAN_deg    = rad2deg(O_);
    H.argp_deg    = rad2deg(w_);

    H.RAAN_un_deg = rad2deg(unwrap(O_));
    H.argp_un_deg = rad2deg(unwrap(w_));

    H.d_RAAN_deg  = H.RAAN_un_deg(end) - H.RAAN_un_deg(1);
    H.d_argp_deg  = H.argp_un_deg(end) - H.argp_un_deg(1);

    H.equatorial  = eqf;
end

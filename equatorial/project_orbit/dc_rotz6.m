function R6 = dc_rotz6(ang)
%DC_ROTZ6  Rotate a full 6x1 state (position AND velocity) about the z-axis.
%
%   R6 = DC_ROTZ6(ang)      ang in radians, counter-clockwise
%
%   R6 = [Rz 0; 0 Rz], so R6*X rotates both r and v by the same angle.
%   Used to compare X(T_r) with X0 rotated by the apsidal advance, and to
%   build the monodromy matrix in the frame that turns with the ellipse.

c  = cos(ang);  s = sin(ang);
Rz = [c -s 0; s c 0; 0 0 1];
R6 = blkdiag(Rz, Rz);
end

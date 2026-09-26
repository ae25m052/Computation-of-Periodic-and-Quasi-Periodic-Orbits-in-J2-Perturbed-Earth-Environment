function tag = inc_tag(inc_deg)
%INC_TAG  Filename tag for an inclination: 000, 045, 063p43, 116p57, 180.
%
%   tag = inc_tag(63.434949)     ->   '063p43'
%
%   WHY ZERO PADDING
%   Windows Explorer and MATLAB's dir() both sort file names
%   ALPHABETICALLY, not numerically. Without padding, "inc15" sorts before
%   "inc5" and a folder of sweep figures is impossible to read in order.
%   Three digits fixes that for everything up to 180.
%
%   WHY "p" INSTEAD OF A DECIMAL POINT
%   It keeps the file name to a single dot -- the one before the extension
%   -- and it still sorts correctly: 063p43 falls neatly between 060 and
%   075. (A dot inside a FOLDER name is a separate and worse problem in
%   MATLAB, which reads it as a package separator; this avoids the habit
%   entirely.)
%
%   M.Tech. project, Raj Khismatrao (AE25M052).

    if abs(inc_deg - round(inc_deg)) < 1e-6
        tag = sprintf('%03d', round(inc_deg));
    else
        whole = floor(inc_deg);
        frac  = round((inc_deg - whole) * 100);
        if frac >= 100          % guard the carry, e.g. 62.999 -> 063
            whole = whole + 1;
            frac  = 0;
        end
        tag = sprintf('%03dp%02d', whole, frac);
    end
end

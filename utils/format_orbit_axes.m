function format_orbit_axes(X, ttl)
%FORMAT_ORBIT_AXES  Common axis dressing for every 3D orbit figure.
%
%   format_orbit_axes(X, 'Kepler (\epsilon = 0),  i = 45^\circ,  400 revs')
%
%   Inputs:
%       X     3xN or 6xN state history, used only to size the box
%       ttl   title string (may contain TeX, may be a cell array of lines)
%
%   EQUAL CUBIC LIMITS, DRIVEN BY max|r|
%   This matters more than it looks. MATLAB's autoscaling would fit the box
%   tightly to the data, which squashes the Earth sphere into an ellipsoid
%   and, worse, puts two figures of the same orbit at different
%   inclinations on different scales -- so they cannot be compared side by
%   side. Driving the limits from max|r| keeps every figure in the sweep on
%   one common scale.
%
%   WHY sgtitle AND NOT title
%   A plain axes title sits immediately above the axes box. On a 3D plot
%   with `axis vis3d` the axes box reaches almost to the top of the figure,
%   so once export_for_slides enlarges the fonts the title overflows the
%   figure and exportgraphics crops its top row of pixels off. sgtitle is a
%   FIGURE-level title: MATLAB reserves space for it by shrinking the axes,
%   so it cannot be clipped no matter how large the font becomes.
%
%   Operates on the current axes and figure, so call it after the orbit has
%   been plotted (and after any colorbar) and before style_light.
%
%   M.Tech. project, Raj Khismatrao (AE25M052).

    grid on;
    set(gca, 'GridAlpha', 0.25);
    view(45, 25);
    axis vis3d;

    xlabel('x  [DU]');
    ylabel('y  [DU]');
    zlabel('z  [DU]');

    lim = 1.10 * max(vecnorm(X(1:3,:)));
    xlim([-lim lim]);
    ylim([-lim lim]);
    zlim([-lim lim]);

    sgtitle(ttl);
end

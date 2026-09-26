function fig = dc_light_figure(pos)
%DC_LIGHT_FIGURE  New figure that stays white even when MATLAB is in dark mode.
%
%   fig = DC_LIGHT_FIGURE([left bottom width height])
%
%   Newer MATLAB versions (R2025a on) have a dark theme. In dark theme a
%   plain  figure('Color','w')  still gets dark axes, so black curves and
%   grey text become invisible. Two fixes, applied together:
%     1. If this MATLAB has the figure 'Theme' property, set it to 'light'.
%        (Older versions do not have it; isprop then returns false and the
%        line is skipped, so the file works on any version.)
%     2. Every axes is also coloured explicitly by dc_style_axes.m.
%
%   'InvertHardcopy','off' makes the saved PNG keep exactly these colours.

fig = figure('Position', pos);
if isprop(fig, 'Theme')
    fig.Theme = 'light';
end
set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
end

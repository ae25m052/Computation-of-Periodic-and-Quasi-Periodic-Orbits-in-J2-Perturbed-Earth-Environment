function eq_style_axes(ax)
%EQ_STYLE_AXES  White background, dark text, light grid, for one axes.
%
%   EQ_STYLE_AXES(ax)
%
%   Call it AFTER plotting, labelling and adding the legend, because it
%   also recolours the title, axis labels, any text() labels and the legend.
%   Colours are set explicitly so the figure looks the same in light and
%   dark MATLAB themes and in the saved PNG.

ink  = [0.10 0.10 0.10];     % text and axis lines
grid_col = [0.82 0.82 0.82];

set(ax, 'Color', 'w', 'XColor', ink, 'YColor', ink, ...
        'GridColor', grid_col, 'GridAlpha', 1, 'LineWidth', 0.8, ...
        'FontSize', 11, 'Box', 'on');
grid(ax, 'on');
ax.Title.Color  = ink;
ax.XLabel.Color = ink;
ax.YLabel.Color = ink;

txt = findobj(ax, 'Type', 'text');
set(txt, 'Color', ink);

lg = ax.Legend;                      % empty if this axes has no legend
if ~isempty(lg)
    set(lg, 'Color', 'w', 'TextColor', ink, 'EdgeColor', grid_col);
end
end

function export_for_slides(figHandle, filename, opts)
%EXPORT_FOR_SLIDES Restyle a figure for projection and export it at high DPI.
%
%   export_for_slides(gcf, 'myfig.png')
%   export_for_slides(gcf, 'myfig.png', struct('FontSize',20,'DPI',300))
%
%   WHY THIS EXISTS
%   The review template requires that "all the letters/numbers in the
%   figures/tables used should be similar in size to the font size" of the
%   body text -- i.e. >= 20pt Calibri equivalent. MATLAB's default axes
%   font is ~10pt, so a figure exported straight from the screen and then
%   shrunk onto a slide has axis labels roughly a quarter of the required
%   size. The template devotes an entire slide ("A BAD FIGURE") to this
%   exact failure, so it is worth fixing at the source.
%
%   It also forces a light background with DARK text. If MATLAB is running
%   in dark theme, axes text defaults to near-white; exporting that onto a
%   white slide gives pale-grey-on-white labels that vanish when projected.
%
%   Inputs:
%       figHandle - figure handle, e.g. gcf
%       filename  - output file, e.g. 'kep_vs_j2.png'
%       opts      - optional struct:
%                     .FontSize   base font size in points   (default 20)
%                     .TitleSize  title font size            (default 22)
%                     .LineWidth  trace line width           (default 2.0)
%                     .DPI        export resolution          (default 300)
%                     .Width      exported width in inches   (default 10)
%                     .Height     exported height in inches  (default 6)
%
%   NOTE: this modifies the figure that is passed in. Re-run the
%   generating script if you want the original screen styling back.

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'FontSize'),  opts.FontSize  = 20;  end
    if ~isfield(opts, 'TitleSize'), opts.TitleSize = 22;  end
    if ~isfield(opts, 'LineWidth'), opts.LineWidth = 2.0; end
    if ~isfield(opts, 'DPI'),       opts.DPI       = 300; end
    if ~isfield(opts, 'Width'),     opts.Width     = 10;  end
    if ~isfield(opts, 'Height'),    opts.Height    = 6;   end

    dark = [0.10 0.10 0.10];

    set(figHandle, 'Color', 'w', 'InvertHardcopy', 'off');

    ax = findall(figHandle, 'Type', 'axes');
    for k = 1:numel(ax)
        a = ax(k);
        set(a, 'FontSize',  opts.FontSize, ...
               'FontName',  'Calibri', ...
               'LineWidth', 1.2, ...
               'XColor', dark, 'YColor', dark, 'ZColor', dark, ...
               'GridColor', [0.35 0.35 0.35], 'GridAlpha', 0.30);

        % Axis labels and title: dark, and a touch larger than tick labels
        set([a.XLabel a.YLabel a.ZLabel], ...
            'FontSize', opts.FontSize, 'Color', dark, 'FontName', 'Calibri');
        set(a.Title, 'FontSize', opts.TitleSize, 'Color', dark, ...
            'FontName', 'Calibri', 'FontWeight', 'bold');

        % Thicken plotted lines so they survive projection
        ln = findall(a, 'Type', 'line');
        for j = 1:numel(ln)
            if get(ln(j), 'LineWidth') < opts.LineWidth
                set(ln(j), 'LineWidth', opts.LineWidth);
            end
        end
    end

    % Legends and colorbars carry their own font settings
    lg = findall(figHandle, 'Type', 'legend');
    set(lg, 'FontSize', max(opts.FontSize - 2, 14), 'FontName', 'Calibri', ...
            'TextColor', dark);
    cb = findall(figHandle, 'Type', 'colorbar');
    set(cb, 'FontSize', max(opts.FontSize - 2, 14), 'Color', dark);

    % Overall figure titles (sgtitle) are annotation-like text objects
    tx = findall(figHandle, 'Type', 'text');
    for k = 1:numel(tx)
        if isempty(get(tx(k), 'Parent')), continue; end
        c = get(tx(k), 'Color');
        if all(c > 0.6)               % was light-on-dark; flip to dark
            set(tx(k), 'Color', dark);
        end
    end

    % Fixed physical size so exported text size is predictable
    set(figHandle, 'Units', 'inches', ...
                   'Position', [1 1 opts.Width opts.Height], ...
                   'PaperUnits', 'inches', ...
                   'PaperPosition', [0 0 opts.Width opts.Height], ...
                   'PaperSize', [opts.Width opts.Height]);

    drawnow;
    exportgraphics(figHandle, filename, 'Resolution', opts.DPI, ...
                   'BackgroundColor', 'white');
    fprintf('Exported %s at %d DPI (%.1f x %.1f in, base font %dpt)\n', ...
            filename, opts.DPI, opts.Width, opts.Height, opts.FontSize);
end

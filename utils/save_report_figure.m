function save_report_figure(fig, filename, width_cm, height_cm, fontsize_pt)
%SAVE_REPORT_FIGURE  Export a figure at a FIXED physical size and font size.
%
%   save_report_figure(gcf, 'fig.png', width_cm, height_cm)
%   save_report_figure(gcf, 'fig.png', width_cm, height_cm, fontsize_pt)
%
%   WHY THIS EXISTS
%   Text inside a figure is rendered at whatever size MATLAB used, then
%   LaTeX rescales the whole image to fit the page. If two figures are
%   exported at different sizes and then placed at the same size on the
%   page, they are rescaled by different factors and their fonts come out
%   different -- even if the MATLAB font size was identical.
%
%   The fix is to export at the size you will actually use, and then in
%   LaTeX include it at THAT SAME SIZE. Nothing gets rescaled, so the
%   fonts land exactly as specified.
%
%   THE SIZES THIS REPORT EXPECTS (A4, 1 inch margins -> 15.9 cm of text)
%
%     Figure 1  Keplerian_vs_J2-Perturbed_Orbit   15.9 x  8.6 cm
%     Figure 2  Keplerian_trajectory               6.0 x  6.0 cm
%     Figure 2  J2_trajectory_coloured_per_orbit   8.0 x  6.0 cm
%     Figure 3  J2_minus_Kepler                   14.0 x  7.6 cm
%
%   fontsize_pt defaults to 9. The report body is 11 pt; figure text one
%   or two points smaller reads correctly and is the usual convention.
%
%   EXAMPLE
%       main_stage1_keplerian_seed        % draws the figure
%       title('')                         % drop the in-figure title,
%                                         % the LaTeX caption says it
%       save_report_figure(gcf, 'Keplerian_trajectory.png', 6.0, 6.0)

if nargin < 5 || isempty(fontsize_pt), fontsize_pt = 9;   end
if nargin < 4 || isempty(height_cm),   height_cm   = 6.0; end
if nargin < 3 || isempty(width_cm),    width_cm    = 8.0; end
if isempty(fig), fig = gcf; end

% ------------------------------------------------------------------
% 1. Set the figure to the exact physical size it will occupy.
% ------------------------------------------------------------------
set(fig, 'Units', 'centimeters');
pos = get(fig, 'Position');
set(fig, 'Position', [pos(1) pos(2) width_cm height_cm]);

set(fig, 'PaperUnits', 'centimeters', ...
         'PaperPosition', [0 0 width_cm height_cm], ...
         'PaperSize', [width_cm height_cm]);

% ------------------------------------------------------------------
% 2. Force ONE font size and family on every text object in the figure:
%    axes tick labels, axis labels, title, legend, colorbar, annotations.
%    findall (not findobj) is needed because some of these are hidden
%    handles.
% ------------------------------------------------------------------
set(findall(fig, '-property', 'FontSize'),   'FontSize',   fontsize_pt);
set(findall(fig, '-property', 'FontName'),   'FontName',   'Helvetica');
set(findall(fig, '-property', 'FontWeight'), 'FontWeight', 'normal');

% Colorbar labels sometimes carry their own multiplier; reset it.
cbs = findall(fig, 'Type', 'ColorBar');
for k = 1:numel(cbs)
    cbs(k).FontSize = fontsize_pt;
    cbs(k).Label.FontSize = fontsize_pt;
    cbs(k).Label.FontWeight = 'normal';
end

% ------------------------------------------------------------------
% 3. Export. exportgraphics crops surrounding whitespace automatically,
%    which saveas and print do not.
% ------------------------------------------------------------------
exportgraphics(fig, filename, ...
    'Resolution', 300, ...
    'BackgroundColor', 'white');

% ------------------------------------------------------------------
% 4. Report what was produced, so you can check it against the LaTeX.
% ------------------------------------------------------------------
info = imfinfo(filename);
fprintf('%s\n', filename);
fprintf('  requested : %.2f x %.2f cm at %d pt\n', width_cm, height_cm, fontsize_pt);
fprintf('  written   : %d x %d px  (%.0f px/cm)\n', ...
        info.Width, info.Height, info.Width/width_cm);
fprintf('  LaTeX     : \\includegraphics[width=%.1fcm]{%s}\n\n', ...
        width_cm, filename);

end

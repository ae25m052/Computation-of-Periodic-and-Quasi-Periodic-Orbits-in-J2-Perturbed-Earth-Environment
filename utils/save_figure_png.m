function save_figure_png(figh, folder, name, fopts, figfolder)
%SAVE_FIGURE_PNG  Export one figure as PNG, and optionally as a MATLAB .fig.
%
%   save_figure_png(fig, dirJ, 'j2_inc045', struct('Width',10,'Height',8))
%   save_figure_png(fig, dirJ, 'j2_inc045', opts, dirFJ)   % also writes .fig
%
%   Inputs:
%       figh      figure handle
%       folder    destination folder for the PNG (created if it does not exist)
%       name      file name WITHOUT the extension
%       fopts     optional struct passed straight to export_for_slides
%                 (.FontSize .TitleSize .LineWidth .DPI .Width .Height)
%       figfolder optional folder for a MATLAB .fig copy of the same figure.
%                 Leave it out, or pass '', to write only the PNG.
%
%   The PNG goes through the project's export_for_slides.m, which enforces
%   the review template's >= 20 pt requirement on every label and exports at
%   300 DPI on a white background with dark text.
%
%   WHY A .fig AS WELL
%   A PNG is a flat picture. A .fig keeps the actual graphics objects, so a
%   3D orbit saved this way can be reopened, ROTATED, zoomed, recoloured,
%   relabelled and re-exported -- none of which is possible from the PNG.
%   The .fig is written AFTER the PNG export, so the two match exactly.
%
%   THE VISIBILITY GOTCHA
%   savefig stores the figure's Visible property, and openfig honours it. A
%   figure built off-screen (SHOW_FIGURES = false, which is what a long
%   sweep wants) would therefore be saved as invisible, and double-clicking
%   the .fig later would appear to do nothing at all. The Visible flag is
%   set back on immediately before saving so the file reopens normally.
%
%   Folders are created here rather than only at the start of the run, so a
%   missing or mistyped folder cannot destroy a sweep that has already spent
%   twenty minutes integrating.
%
%   NOTE: export_for_slides MODIFIES the figure it is handed (fonts, colours
%   and physical size). So this should be the LAST thing done to a figure
%   before it is closed.
%
%   M.Tech. project, Raj Khismatrao (AE25M052).

    if nargin < 4 || isempty(fopts)
        fopts = struct();
    end
    if nargin < 5
        figfolder = '';
    end

    % ---- PNG ----------------------------------------------------------
    if ~exist(folder, 'dir')
        [ok, msg] = mkdir(folder);
        if ~ok
            error('save_figure_png:mkdir', ...
                  'Could not create folder %s: %s', folder, msg);
        end
    end

    export_for_slides(figh, fullfile(folder, [name '.png']), fopts);

    % ---- MATLAB .fig --------------------------------------------------
    if ~isempty(figfolder)
        if ~exist(figfolder, 'dir')
            [ok, msg] = mkdir(figfolder);
            if ~ok
                error('save_figure_png:mkdir', ...
                      'Could not create folder %s: %s', figfolder, msg);
            end
        end
        set(figh, 'Visible', 'on');      % see THE VISIBILITY GOTCHA above
        savefig(figh, fullfile(figfolder, [name '.fig']));
    end
end

function sol = extend_uvabo_ransac(sol, varargin)
% EXTEND_UVABO_RANSAC Extend a relaxed TDOA solution
%   [solout, inlier_count] = EXTEND_UVABO_RANSAC(sol) extends a relaxed
%       solution to utilize more rows and columns of the TDOA measurements.
%   [solout, inlier_count] = EXTEND_UVABO_RANSAC(sol, ...) additional
%       inputs:
%           display - controls the amount of printing.
%               off, none - no printing.
%               iter - printouts for each iteration of RANSAC/optimization.
%           iters - number of iterations in RANSAC loop.
%           min_inliers - minimum number of inliers in new row or column.
%           refine_iters - number of iterations of local optimization
%               between each extension.
%           threshold - threshold for the absolute error in TDOA measurment
%               when classifying inliers/outliers.

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'iters', 500);
addParameter(p, 'min_inliers', 3);
addParameter(p, 'refine_iters', 5);
addParameter(p, 'threshold', 0.2);
parse(p, varargin{:});
opts = p.Results;

ext_opts.iters = opts.iters;
ext_opts.min_inliers = opts.min_inliers;
ext_opts.threshold = opts.threshold;

extended = 1;
while extended
    [sol1, count1] = extend_ua_ransac(sol, ext_opts);
    [sol2, count2] = extend_vbo_ransac(sol, ext_opts);

    if count1 == 0 && count2 == 0
        extended = 0;
    else
        if count1 > count2
            sol = sol1;
        else
            sol = sol2;
        end

        sol = refine_uvabo(sol, 'max_iters', opts.refine_iters, ...
            'display', opts.display);

        if strcmpi(opts.display, 'iter') % TODO: Option to enable plots?
            misstdoa_briefer_report(sol);
            misstdoa_brief_visualization(sol);
        end
    end
end

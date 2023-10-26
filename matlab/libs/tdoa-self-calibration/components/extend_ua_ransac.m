function [solout, inlier_count] = extend_ua_ransac(sol, varargin)
% EXTEND_UA_RANSAC Extend a relaxed TDOA solution over receivers
%   [solout, inlier_count] = EXTEND_UA_RANSAC(sol) extends a relaxed
%       solution with one row in u and a using robust methods.
%   [solout, inlier_count] = EXTEND_UA_RANSAC(sol, ...) additional inputs:
%           iters - number of iterations in RANSAC loop.
%           min_inliers - minimum number of inliers in new row.
%           threshold - threshold for the absolute error in TDOA measurment
%               when classifying inliers/outliers.

% Parse inputs.
p = inputParser;
addParameter(p, 'iters', 100);
addParameter(p, 'min_inliers', 3);
addParameter(p, 'threshold', 0.2);
parse(p, varargin{:});
opts = p.Results;

solout = sol;
inlier_count = 0;

z = sol.z;
max_inliers = 0;
bestsol = sol;

possible_rows = 1:size(z, 1);
possible_rows(sol.rows) = [];

% We need 4 columns for finding unew and anew, and at least 1 for testing.
possible_rows(sum(isfinite(z(possible_rows, sol.cols)), 2) < 5) = [];

if isempty(possible_rows)
    return;
end

for i = 1:opts.iters
    newrow = possible_rows(randi(length(possible_rows)));
    possible_cols = sol.cols(isfinite(z(newrow, sol.cols)));
    use_cols = possible_cols(randperm(length(possible_cols), 4));
    test_cols = setdiff(possible_cols, use_cols);

    % Get indices in sol.cols for use_cols and test_cols.
    colsorder = zeros(1, length(use_cols));
    for k = 1:length(use_cols)
        colsorder(k) = find(sol.cols == use_cols(k));
    end
    colsorder_test = zeros(1, length(test_cols));
    for k = 1:length(test_cols)
        colsorder_test(k) = find(sol.cols == test_cols(k));
    end

    % Estimate new row in u and a using z(newrow, use_cols).
    vsub = sol.v(:, colsorder);
    bsub = sol.b(colsorder);
    osub = sol.o(colsorder);
    zsub = z(newrow, use_cols);
    AA = (zsub - osub).^2 - bsub;
    bb = [-2 * vsub; ones(1, length(use_cols))];
    xx = AA / bb;
    unew = xx(1:3);
    anew = xx(4);

    % Evaluate using columns in test_cols.
    vtest = sol.v(:, colsorder_test);
    btest = sol.b(colsorder_test);
    otest = sol.o(colsorder_test);
    err = z(newrow, test_cols) - (sqrt(-2 * unew * vtest + anew + btest) + otest);

    inind = abs(err) < opts.threshold;
    nrin = sum(inind);

    if nrin > max_inliers
        max_inliers = nrin;
        bestsol = sol;
        bestsol.rows = [bestsol.rows, newrow];
        bestsol.u = [bestsol.u; unew];
        bestsol.a = [bestsol.a; anew];
        inlcols = [use_cols, test_cols(inind)];
        bestsol.inlmatrix(newrow, inlcols) = 1;
    end
end

if max_inliers >= opts.min_inliers
    solout = bestsol;
    inlier_count = max_inliers;
end


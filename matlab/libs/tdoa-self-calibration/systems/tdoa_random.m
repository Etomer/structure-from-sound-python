function [r, s, o, sol] = tdoa_random(z, varargin)
% TDOA Perform TDOA self-calibration using a naive method
%   [r, s, o, sol] = TDOA(z) performs TDOA self-calibration in 3D using the
%       TDOA measurements in z (m x n). The estimated receiver (3 x m) and
%       sender (3 x n) positions, and the offsets (1 x n) are returned
%       together with a solution structure containing additional
%       information. The method works by randomly initialize the node
%       positions, followed by local optimization.
%   [r, s, o, sol] = TDOA(z, ...) additional inputs:
%       display - controls the amount of printing.
%           off, none - no printing.
%           iter - printouts for each iteration of RANSAC/optimization.
%       sigma - estimated standard deviation of the measurement noise. This
%           is used to set suitable thresholds for RANSAC.
%       offsetsolver - a structure with the offset solver to use.
%
% See also get_offset_solvers.

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'sigma', 0.01);
addParameter(p, 'inits', 1);
addParameter(p, 'mul5', inf);
parse(p, varargin{:});
opts = p.Results;

% TODO: Generalize for other dimensions than 3D.

besterr = inf;
for i=1:opts.inits
    % Find initial r, s, o.
    sol = init_rso_random(z);

    % Local optimization over r, s, o.
    sol = refine_rso(sol, 'display', opts.display);
    if isfinite(opts.mul5)
        sol = refine_rso_robust(sol, 'display', opts.display, 'threshold', opts.mul5*opts.sigma);
    end

    % Calculate residuals.
    d = pdist2(sol.r', sol.s');
    res = sol.z(sol.rows, sol.cols) - sol.o - d;
    inlres = res(sol.inlmatrix(sol.rows, sol.cols));
    err = rms(inlres);
    if err < besterr
        besterr = err;
        bestsol = sol;
    end
end
sol = bestsol;

% Output.
[m, n] = size(z);
dim = 3;
r = nan(dim, m);
s = nan(dim, n);
o = nan(1, n);

r(:, sol.rows) = sol.r;
s(:, sol.cols) = sol.s;
o(:, sol.cols) = sol.o;

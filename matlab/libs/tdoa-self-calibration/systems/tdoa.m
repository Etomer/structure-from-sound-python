function [r, s, o, sol] = tdoa(z, varargin)
% TDOA Perform TDOA self-calibration using robust methods
%   [r, s, o, sol] = TDOA(z) performs TDOA self-calibration in 3D using the
%       TDOA measurements in z (m x n). The estimated receiver (3 x m) and
%       sender (3 x n) positions, and the offsets (1 x n) are returned
%       together with a solution structure containing additional
%       information. The method works as follows:
%           1. Initialize a relaxed solution in u, v, a, b, o.
%           2. Extend solution to more rows and columns.
%           3. Upgrade to solution in r, s, o.
%   [r, s, o, sol] = TDOA(z, ...) additional inputs:
%       display - controls the amount of printing.
%           off, none - no printing.
%           iter - printouts for each iteration of RANSAC/optimization.
%       sigma - estimated standard deviation of the measurement noise. This
%           is used to set suitable thresholds for RANSAC.
%       offsetsolver - a structure with the offset solver to use.
%
% See M. Larsson et al. (2021) Fast and Robust Stratified Self-Calibration
% Using Time-Difference-of-Arrival Measurements for details.
%
% See also get_offset_solvers.

% Use (9r/5s) offset solver by default.
default_solver = get_offset_solvers(3, 9, 5);

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'offsetsolver', default_solver);
addParameter(p, 'sigma', 0.01);
addParameter(p, 'mul1', 4); % Sets threshold based on sigma.
addParameter(p, 'mul2', 6); % Sets threshold based on sigma.
addParameter(p, 'mul3', 10); % Sets threshold based on sigma.
addParameter(p, 'mul4', 12); % Sets threshold based on sigma.
addParameter(p, 'mul5', inf); % Sets threshold based on sigma.
parse(p, varargin{:});
opts = p.Results;

% TODO: Generalize for other dimensions than 3D.

% Find initial u, v, a, b, o.
sol = init_uvabo_ransac(z, ...
    'display', opts.display, 'threshold', opts.mul1*opts.sigma,...
    'solver', opts.offsetsolver);

% Bundle over u, v, a, b, o.
sol = refine_uvabo(sol, 'display', opts.display);

% Extend solutions to more rows and columns.
sol = extend_uvabo_ransac(sol, ...
    'display', opts.display, 'threshold', opts.mul2*opts.sigma);

% Upgrade to solution in r, s, o.
sol = upgrade_ransac(sol, ...
    'display', opts.display, 'threshold', opts.mul3*opts.sigma);
% sol = upgrade_linear(sol);

% Extend solutions to more rows and columns.
% TODO: Clean up and rename functions below.
sol = misstdoa_reestimate_cols_rso(sol, 'threshold', opts.mul4*opts.sigma);
sol = misstdoa_reestimate_rows_rso(sol, 'threshold', opts.mul4*opts.sigma);

% Local optimization over L, q.
% sol = refine_Lq(sol);

% Local optimization over r, s, o.
sol = refine_rso(sol, 'display', opts.display, 'max_iters', 200);
if isfinite(opts.mul5)
    sol = refine_rso_robust(sol, 'display', opts.display, 'threshold', opts.mul5*opts.sigma);
end

% Output.
[m, n] = size(z);
dim = 3;
r = nan(dim, m);
s = nan(dim, n);
o = nan(1, n);

r(:, sol.rows) = sol.r;
s(:, sol.cols) = sol.s;
o(:, sol.cols) = sol.o;

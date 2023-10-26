function solout = upgrade_ransac(sol, varargin)
% UPGRADE_RANSAC Upgrade a relaxed solution
%   sol = UPGRADE_RANSAC(sol) upgrades a relaxed solution, a solution in u,
%       v, a, b and o, to a solution in r, s and o. This is done using
%       robust methods (see note below).
%   sol = UPGRADE_RANSAC(sol, ...) additional inputs:
%       display - controls the amount of printing.
%           off, none - no printing.
%           iter - printouts for each iteration of RANSAC/optimization.
%       iters - number of iterations in RANSAC loop.
%       threshold - threshold for the absolute error in TDOA measurment
%           when classifying inliers/outliers.
%
% See M. Larsson et al. (2020) Upgrade Methods for Stratified Sensor
% Network Self-calibration for details.

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'iters', 50);
addParameter(p, 'threshold', 0.1);
parse(p, varargin{:});
opts = p.Results;

% Display.
if ~any(strcmpi(opts.display, {'off', 'none'}))
    fprintf('Upgrading solution using RANSAC.\n');
end

% Initialize output struct.
solout.type = 'rso';
solout.rows = sol.rows;
solout.cols = sol.cols;
solout.inlmatrix = sol.inlmatrix;
solout.z = sol.z;
solout.o = sol.o;

% Convert sol to the problem structure used in the previous paper.
prob = struct();
prob.a = sol.b;
prob.b = sol.a;
prob.c = 0;
prob.U = sol.u';
prob.V = sol.v;
prob.Dmeas = sol.z(sol.rows, sol.cols) - sol.o;
inliers = sol.inlmatrix(sol.rows, sol.cols);

% Options controlling local optimization in minimal solvers.
solvopts = struct();
solvopts.tol = 1e-9;
solvopts.maxIters = 20;
solvopts.refine = false;

% Get all upgrade solvers.
solvers = getSolvers();

% Keep only those who work with the number of receivers and senders.
allids = vertcat(solvers.id);
keep = length(sol.rows) >= allids(:, 1) + 1 & length(sol.cols) >= allids(:, 2) + 1;
solvers = solvers(keep);

% RANSAC loop.
max_inliers = 0;
for i = 1:opts.iters
    % Pick a solver randomly.
    % TODO: Is this the best strategy?
    solver = solvers(randi(length(solvers)));

    % Create random minimal sample from the problem. Since u and v are
    % dense we do not need to worry about missing data.
    sample = createRandomSample(solver.id, prob);

    % Solve minimal problem.
    [Lti, q] = solver.solve(sample, solvopts);

    for j = 1:length(Lti)
        Rhat = Lti{j} * sample.fullU;
        Shat = Lti{j}' \ (sample.fullV + q(:, j));

        % Calculate error in distances.
        Dhat = pdist2(Rhat', Shat');
        err = abs(prob.Dmeas(inliers)-Dhat(inliers));
        ninl = sum(err < opts.threshold);

        if ninl > max_inliers
            max_inliers = ninl;
            solout.r = Rhat;
            solout.s = Shat;
            solout.L = inv(Lti{j}');
            solout.q = q(:, j);
            solout.u = sample.fullU';
            solout.v = sample.fullV;

            if strcmpi(opts.display, 'iter')
                fprintf('Iter %3d: inliers = %3d, solver = %s\n',...
                    i, ninl, solver.name);
            end
        end
    end
end

if max_inliers == 0
    error('RANSAC upgrade failed to find a solution.');
end

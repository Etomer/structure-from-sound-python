function [bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, varargin)
% INIT_UVABO_RANSAC Initialize a relaxed TDOA solution
%   sol = INIT_UVABO_RANSAC(z) initializes a relaxed solution, i.e., a
%       solution in u, v, a, b and o, using robust methods.
%   sol = INIT_UVABO_RANSAC(z) additional inputs:
%           display - controls the amount of printing.
%               off, none - no printing.
%               iter - printouts for each iteration of RANSAC/optimization.
%           iters - number of iterations in RANSAC loop.
%           solver - an array of offset solvers to use.
%           threshold - threshold for the absolute error in TDOA measurment
%               when classifying inliers/outliers.

% 1. The offsets are solved for first using the minimal offset solvers.
% 2. The minimal solution is then extended to include more columns using
%    robust methods.
% 3. Inliers are counted and the best solution is kept.

% Default solver.
default_solver.solv = @solver_toda_rank3_95;
default_solver.name = func2str(default_solver.solv);
default_solver.m = 9;
default_solver.n = 5;
default_solver.rank = 3;

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'iters', 100000);
addParameter(p, 'solver', default_solver);
%addParameter(p, 'threshold', 0.2);
addParameter(p, 'threshold', 0.01);
parse(p, varargin{:});
opts = p.Results;

% Display.
if ~any(strcmpi(opts.display, {'off', 'none'}))
    fprintf('Running RANSAC over U, V, a, b, o.\n');
end

% Init RANSAC.
bestsol = [];
max_inliers = 0;
best_err = Inf;
stats1 = zeros(1, opts.iters);
stats2 = zeros(3, 0);
stats2_counter = 0;

for iRansac = 1:opts.iters
    % m receivers and n senders.
    [m, n] = size(z);
    ok = isfinite(z);
    
    % Pick random solver.
    solver = opts.solver(randi(length(opts.solver)));

    % m_solv, n_solv are the number of rows and columns needed for the minimal
    % solver.
    m_solv = solver.m;
    n_solv = solver.n;
    rank_solv = solver.rank;

    % Precompute matrices needed for going from (z,o) to (u,v,a,b);
    wr = zeros(m_solv, 1);
    wr(1) = 1;
    ws = zeros(n_solv, 1);
    ws(1) = 1;

    Cr = eye(m_solv) - wr * ones(1, m_solv);
    Cs = eye(n_solv) - ws * ones(1, n_solv);
    
    % Choose m_solv random rows.
    rowsol = randperm(m, m_solv);

    % Select all columns which has data for all rows.
    ok_cols = find(all(ok(rowsol, :), 1));
    if length(ok_cols) < n_solv
        continue;
    end
    colsol = ok_cols(randperm(length(ok_cols), n_solv));

    sols = solver.solv(z(rowsol, colsol));
    nsols = size(sols, 2);
    stats1(1, iRansac) = nsols;
    for iSol = 1:nsols
        osol = sols(:, iSol).';
        if any(abs(imag(osol))./abs(osol) > 1e-6)
            continue;
        end
        osol = real(osol);
        zsol = z(rowsol, colsol);
        dsol = zsol - osol;
        if any(dsol(:) <= 0)
            continue;
        end

        % TODO: Local optimization over offsets?

        % Find initial estimate of u,v,a,b from zsol.
        dsol2 = dsol.^2;
        M = Cr' * dsol2 * Cs / (-2); % This M is M/(-2) compared to the ICASSP 2020 paper
        [uu, ss, vv] = svd(M);
        u = uu(:, 1:rank_solv);
        v = ss(1:rank_solv, 1:rank_solv) * vv(:, 1:rank_solv)';
        a = dsol2 * ws - wr' * dsol2 * ws;
        b = wr' * dsol2;

        % The following should now hold:
        % dsol.^2 = (zsol-osol).^2 = -2*(u*v)+a+b


        % TODO: Clean up code below and extract function.

        % Now check for inliers among the rest of the cols??
        % Which should we try?
        % all cols that are
        % (i) not in pp and
        % (ii) for which there are at least 6 measurements
        % except for those with indices in pp
        okcol = find(sum(ok(rowsol, :)) >= 6);
        restcols = setdiff(okcol, colsol);
        inliersrest = zeros(size(restcols));
        v_rest = zeros(rank_solv, length(restcols));
        b_rest = zeros(1, length(restcols));
        o_rest = zeros(1, length(restcols));
        inl_rest = zeros(m_solv, length(restcols));
        nr_inliers = 0;
        tot_err = 0;
        for ji = 1:length(restcols)
            j = restcols(ji);
            % For each new column,
            % use the fact that we know (u,a)
            % and
            % calculate (vny,bny,ony)
            % so that (z(cc,j)-ony).^2  is equal to -2*(u*vny)+a+bny
            okrow = find(ok(rowsol, j));
            % For a new column we select a minimum number of 5 rows
            sel5 = okrow(randperm(length(okrow), 5));
            z_cut = z(rowsol(sel5), j);
            u_cut = u(sel5, :);
            a_cut = a(sel5);
            AAA = [-2 * z_cut, ones(5, 1), (-u_cut), -ones(5, 1)];
            bbb = a_cut - z_cut.^2;
            x_part = AAA \ bbb;
            x_hom = [0, 1, 0, 0, 0, 1]';
            ony = x_part(1);
            lamb = ony^2 - x_part(2);
            xxx = x_part + lamb * x_hom;
            vny = xxx(3:5) / (-2); % TODO: Fix to get the -2 right, i.e M/(-2)
            bny = xxx(6);
            v_rest(:, ji) = vny;
            b_rest(1, ji) = bny;
            o_rest(1, ji) = ony;
            % ... and then check for inliers among measurements in
            % this column
            err = (sqrt(relu(-2 * (u(okrow, :) * vny) + a(okrow) + bny)) + ony) - z(rowsol(okrow), j);
            inlid = find(abs(err) < opts.threshold);
            if length(inlid) > 5
                inliersrest(ji) = 1;
                inl_rest(okrow(inlid), ji) = ones(length(inlid), 1);
                nr_inliers = nr_inliers + length(inlid) - 5;
                tot_err = tot_err + sum(err(inlid).^2); % I am adding the five zeros here, but nevermind.
            end
        end


        tot_err = sqrt(tot_err);
%         nr_inliers = sum(inliersrest) + n_solv;
        nr_inliers = nnz(inl_rest) + m_solv * n_solv;

        stats2_counter = stats2_counter + 1;
        stats2(1, stats2_counter) = iRansac;
        stats2(2, stats2_counter) = nr_inliers;
        stats2(3, stats2_counter) = tot_err;
        %stats2(4,stats2_counter)=otmp;

        if (nr_inliers > max_inliers) || ((nr_inliers == max_inliers) && (tot_err < best_err))
            if strcmpi(opts.display, 'iter')
                fprintf('Iter %5d: inliers = %3d, error = %e\n', iRansac, nr_inliers, tot_err);
            end
            %keyboard;

            max_inliers = nr_inliers;
            best_err = tot_err;

            inliersrest = find(inliersrest);
            bestsol.rows = rowsol;
            bestsol.cols = [colsol, restcols(inliersrest)];
            bestsol.row1 = bestsol.rows(1);
            bestsol.col1 = bestsol.cols(1);

            %keyboard;
            o = [osol, o_rest(inliersrest)];
            v = [v, v_rest(:, inliersrest)];
            b = [b, b_rest(inliersrest)];

            bestsol.inlmatrix = false(m, n);
            bestsol.inlmatrix(bestsol.rows, bestsol.cols) = [true(m_solv, n_solv), inl_rest(:, inliersrest)];
            bestsol.z = z;
            bestsol.o = o;
            bestsol.b = b;
            bestsol.a = a;
            bestsol.u = u;
            bestsol.v = v;
            bestsol.type = 'uvabo';
        end
    end
end

if isempty(bestsol)
    error('Failed to find a solution.');
end

if strcmpi(opts.display, 'iter') % TODO: Option to enable plots?
    misstdoa_briefer_report(bestsol);
    misstdoa_brief_visualization(bestsol);
end

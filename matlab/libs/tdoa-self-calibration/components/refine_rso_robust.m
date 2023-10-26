function sol = refine_rso_robust(sol, varargin)
% REFINE_RSO_ROBUST Perform local optimization of receiver and senders
%   sol = REFINE_RSO_ROBUST(sol) performs robust local optimization over
%       the receiver and sender positions, and the offsets. The
%       optimization uses the truncated L2 norm, i.e., all residuals over a
%       specified threshold will be clamped.
%   sol = REFINE_RSO_ROBUST(sol, ...) additional inputs:
%       display - controls the amount of printing.
%           off, none - no printing
%           iter - printouts for each iteration of RANSAC/optimization
%       max_iters - maximum number of interation in optimization.
%       threshold - all residuals over this threshold will be clamped.
%       tol - tolerance for RMS error.

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'max_iters', 100);
addParameter(p, 'threshold', inf);
addParameter(p, 'tol', 1e-6);
parse(p, varargin{:});
opts = p.Results;

% Display.
if ~any(strcmpi(opts.display, {'off', 'none'}))
    fprintf('Refining solution over R, S, o.\n');
end

% Extract variables for convenience.
z = sol.z(sol.rows, sol.cols);
r = sol.r;
s = sol.s;
o = sol.o;
% inliers = sol.inlmatrix(sol.rows, sol.cols);
nonmissing = isfinite(z);

[I, J] = find(nonmissing);
Z = z(nonmissing);

for i = 1:opts.max_iters
    [res, jac] = calcresandjac(r, s, o, I, J, Z);
    [inliers, res, jac] = truncres(res, jac, opts.threshold);

    if strcmpi(opts.display, 'iter')
        fprintf('Iter %3d: norm(res)=%e, rms(res)=%e\n', i, norm(res), rms(res));
    end

    if rms(res) < opts.tol
        break;
    end

    % Gauss-Newton step.
    dz = -(jac' * jac + 1e-6 * speye(size(jac, 2))) \ (jac' * res);

    [rnew, snew, onew] = update(r, s, o, dz);
    resnew = calcresandjac(rnew, snew, onew, I, J, Z);
    [inliersnew, resnew] = truncres(resnew, [], opts.threshold);
    
    % If no improvement, try reducing the step size.
    j = 0;
    while norm(resnew) >= norm(res)
        dz = dz / 2;
        [rnew, snew, onew] = update(r, s, o, dz);
        resnew = calcresandjac(rnew, snew, onew, I, J, Z);
        [inliersnew, resnew] = truncres(resnew, [], opts.threshold);
        j = j + 1;
        if j > 50
            break;
        end
    end

    if norm(resnew) < norm(res)
        r = rnew;
        s = snew;
        o = onew;
        inliers = inliersnew;
    else
        if strcmpi(opts.display, 'iter')
            fprintf('Stalled\n');
        end
        break;
    end
end

inlmatrix = nonmissing;
inlmatrix(nonmissing) = inliers;

sol.r = r;
sol.s = s;
sol.o = o;
sol.inlmatrix(sol.rows, sol.cols) = inlmatrix;
end

function [inliers, res, jac] = truncres(res, jac, threshold)
% Set threshold to some quantile, i.e., assume outlier ratio.
% [~,sortind] = sort(abs(res(:)),'descend');
% threshold = abs(res(sortind(max(1, round(0.05*length(sortind))))));

resm = res < -threshold;
resp = res > threshold;
inliers = ~resm & ~resp;

% TOOD: Remove debug print.
% fprintf('Truncated residuals: %d (%d %%)\n', nnz(resm | resp), round(100*nnz(resm | resp)/length(res)));

if nargout > 1
    res(resm) = -threshold;
    res(resp) = threshold;
end
if nargout > 2
    jac(~inliers, :) = 0;
end
end

function [res, jac] = calcresandjac(r, s, o, I, J, Z)
v = r(:, I) - s(:, J);
d = vecnorm(v);
res = d' + o(J)' - Z;

if nargout > 1
    dim = size(r, 1);
    m = size(r, 2);
    n = size(s, 2);
    nres = length(Z);

    vid = v ./ d;
    II = repelem((1:nres)', 2*dim+1, 1);
    JJ = [dim * I + (1 - dim:0), dim * m + dim * J + (1 - dim:0), dim * (m + n) + J]';
    VV = [vid; -vid; ones(1, nres)];
    jac = sparse(II(:), JJ(:), VV(:), nres, dim*(m + n)+n);
end
end

function [r, s, o] = update(r, s, o, dz)
dim = size(r, 1);
m = size(r, 2);
n = size(s, 2);

dzr = dz(1:dim*m);
dzs = dz(dim*m+1:dim*(m + n));
dzo = dz(dim*(m + n)+1:dim*(m + n)+n);

r(:) = r(:) + dzr;
s(:) = s(:) + dzs;
o(:) = o(:) + dzo;
end

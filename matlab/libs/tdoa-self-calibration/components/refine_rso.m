function sol = refine_rso(sol, varargin)
% REFINE_RSO Perform local optimization of receiver and senders
%   sol = REFINE_RSO(sol) performs local optimization over the receiver and
%       sender positions, and the offsets.
%   sol = REFINE_RSO(sol, ...) additional inputs:
%       display - controls the amount of printing.
%           off, none - no printing
%           iter - printouts for each iteration of RANSAC/optimization
%       max_iters - maximum number of interations in optimization.
%       tol - tolerance for RMS error.

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'max_iters', 100);
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
inliers = sol.inlmatrix(sol.rows, sol.cols);

[I, J] = find(inliers);
Z = z(inliers);

for i = 1:opts.max_iters
    [res, jac] = calcresandjac(r, s, o, I, J, Z);

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

    % If no improvement, try reducing the step size.
    j = 0;
    while norm(resnew) >= norm(res)
        dz = dz / 2;
        [rnew, snew, onew] = update(r, s, o, dz);
        resnew = calcresandjac(rnew, snew, onew, I, J, Z);
        j = j + 1;
        if j > 50
            break;
        end
    end

    if norm(resnew) < norm(res)
        r = rnew;
        s = snew;
        o = onew;
    else
        if strcmpi(opts.display, 'iter')
            fprintf('Stalled\n');
        end
        break;
    end
end

sol.r = r;
sol.s = s;
sol.o = o;
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

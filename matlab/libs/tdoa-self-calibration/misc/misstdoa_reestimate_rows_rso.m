function [solut] = misstdoa_reestimate_rows_rso(sol, varargin)

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'iters', 10);
addParameter(p, 'threshold', 0.1);
parse(p, varargin{:});
opts = p.Results;

% Display.
if ~any(strcmpi(opts.display, {'off', 'none'}))
    fprintf('Reestimating receiver positions\n');
end

z = sol.z;
s = sol.s;
o = sol.o;
cols = sol.cols;

for ii = 1:size(z,1)
    d = z(ii, cols)' - o';

    [xx, inlny, nr_inliers] = toa_trilateration_one_ransac(d, s, opts.iters, opts.threshold);
    if nr_inliers < 3
        continue;
    end
    
    %[yy,oo,inlny,nr_inlierid,err_rms] = tdoa_trilateration_y_one_ransac(z(rows,jj),r,ransac_k,ransac_tol);
    xx = real(xx);
    [rny, resny] = toa_trilateration_one_bundle(d, s, xx, inlny);
    %    [rny,resny] = tdoa_trilateration_y_one_bundle(z(rows,jj),r,yy,oo,inlny);
    nr_inliers_ny = length(inlny);
    inlny2 = false(1, length(cols));
    inlny2(inlny) = true;
    rms_ny = sqrt(resny'*resny);
    %[length(inlny) std(resny)]

    % Is there a previous guess?
    tmp = find(sol.rows == ii);
    if length(tmp) < 1,
        % No previous estimate of this column
        tmp = size(sol.r, 2) + 1;
        sol.r(:, tmp) = rny;
        sol.rows(tmp) = ii;
        sol.inlmatrix(ii, cols) = inlny2;
    elseif length(tmp) == 1,
        % There is a previous estimate
        r0 = sol.r(:, tmp);
        inl0 = sol.inlmatrix(ii, :);

        zproj = tdoa_calc_u_from_xyo(r0, s, o);
        res0 = (zproj - z(ii, cols))';
        res0 = res0(find(sol.inlmatrix(ii, cols)));
        nr_inliers_0 = sum(inl0);
        rms_0 = sqrt(res0'*res0);
        %[sum(inl0) length(res0) std(res0)]

        if (nr_inliers_ny > nr_inliers_0) || ((nr_inliers_ny == nr_inliers_0) && (rms_ny < rms_0))
            sol.r(:, tmp) = rny;
            sol.inlmatrix(ii, cols) = inlny2';
        end
    else
        % Det h?r vore konstigt.
        error('');
    end

end;

solut = sol;


if 0,
    if nargin < 3,
        y0 = NaN * ones(x_dim, n);
    end

    if nargin < 4,
        o0 = NaN * ones(1, n);
    end

    if nargin < 5,
        index = find(ones(1, n));
    end;

    if nargin < 6,
        inliers = isfinite(u);
    end;

    inl = zeros(size(inliers));

    ransac_tol = 0.02;
end;

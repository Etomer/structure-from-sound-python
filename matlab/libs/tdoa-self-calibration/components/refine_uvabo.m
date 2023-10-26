function [solout, res, jac] = refine_uvabo(sol, varargin)
% REFINE_UVABO Perform local optimization over a relaxed solution
%   sol = REFINE_UVABO(sol) performs local optimization over a relaxed
%       solution.
%   sol = REFINE_UVABO(sol, ...) additional inputs:
%       display - controls the amount of printing.
%           off, none - no printing
%           iter - printouts for each iteration of RANSAC/optimization
%       max_iters - maximum number of interations in optimization.
%       tol - tolerance for norm of residuals.

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'max_iters', 10);
addParameter(p, 'tol', 1e-6); % TODO: Use in code below. RMSE instead?
parse(p, varargin{:});
opts = p.Results;

%keyboard;
if 1,
    jacresfunction = @misstdoa_jacres_v2_fixneg;
else
    jacresfunction = @misstdoa_jacres_v2;
end


% Display.
if ~any(strcmpi(opts.display, {'off', 'none'}))
    fprintf('Refining solution over U, V, a, b, o.\n');
end

% Regularization of Jacobian inverse.
regularization = 1e-2;

% Certain variables are not allowed to change
m1 = length(sol.rows);
n1 = length(sol.cols);
% residual is size Mx1
% parameter order
% w1 w2 w3 v1 v2 v3 c d o
N = 3 * m1 + 3 * n1 + m1 + n1 + n1;
%
Jw1 = 1:m1;
Jw2 = (1:m1) + m1;
Jw3 = (1:m1) + 2 * m1;
Jv1 = (1:n1) + 3 * m1;
Jv2 = (1:n1) + 3 * m1 + 1 * n1;
Jv3 = (1:n1) + 3 * m1 + 2 * n1;
Jc = (1:m1) + 3 * m1 + 3 * n1;
Jd = (1:n1) + 4 * m1 + 3 * n1;
Jo = (1:n1) + 4 * m1 + 4 * n1;
dontmoveindex = [ ...
    Jw1(1:4), Jw2(1:4), Jw3(1:4), ...
    Jv1(1), Jv2(1), Jv3(1), ...
    Jc(1), ...
    ];
EE = speye(N, N);
EE(:, dontmoveindex) = [];

%%
%keyboard;
for kkk = 1:opts.max_iters
    [res, jac] = jacresfunction(sol);

    %dz = -(jac\res);
    %dz = -(jac'*jac+eye(size(jac,2)))\(jac'*res);
    %[u,s,v]=svd(full(jac),0);
    %keyboard;
    jac0 = jac * EE;
    if regularization == 0
        dz = -EE * (jac0 \ res);
    else
        dz = -EE * ((jac0' * jac0 + regularization * speye(size(jac0, 2))) \ (jac0' * res));
    end
    soln = misstdoa_update_v2(sol, dz);
    res2 = jacresfunction(soln);

    %     aa = [norm(res), norm(res + jac * dz), norm(res2)];
    %     bb = aa;
    %     bb = bb - bb(2);
    %     bb = bb / bb(1);
    % Check that the error actually gets smaller
    if norm(res) < norm(res2)
        cc = norm(jac*dz) / norm(res);
        % bad
        % check that the update is sufficiently big
        % otherwise it is maybe just numerical limitations
        if cc > 1e-4
            % OK then it is probably just the linearization that
            % is not good enough for this large step size, decrease
            kkkk = 1;
            while kkkk < 10 && norm(res) < norm(res2)
                dz = dz / 2;
                soln = misstdoa_update_v2(sol, dz);
                res2 = jacresfunction(soln);
                kkkk = kkkk + 1;
            end
        end
    end
    if strcmpi(opts.display, 'iter')
        fprintf('Iter %3d: norm(res) = %e\n', kkk, norm(res2));
    end

    % TODO: What are these prints? Needed? Informative?
    %         aa = [norm(res), norm(res + jac * dz), norm(res2)];
    %         bb = aa;
    %         bb = bb - bb(2);
    %         bb = bb / bb(1);
    %         cc = norm(jac*dz) / norm(res);
    %         %keyboard;
    %         aa
    %         bb
    %         cc
    % Explanation. These can be used while debugging or understanding
    % how the optimization works.
    % aa contains 
    %   1. the norm of the old residuals
    %   2. the norm of the new residuals as estimated using linearization
    %   3. the norm of the new residual
    % Ideally 2 should be lower than 1
    % and 3 should be similar to 2 
    % and therefore also lower than 1
    % If 2 and 3 are almost equal and lower than 1
    % then the linearisation is good and the optimization worked
    % If 2 and 3 are very different then the linearization is bad. 
    % If 3 is higher than 1, then one should make the search step smaller
    % 
    
    if norm(res2) < norm(res)
        sol = soln;
        if norm(res2) < opts.tol
            if strcmpi(opts.display, 'iter')
                fprintf('Residual smaller than tolerance %e. Aborting.\n', opts.tol);
            end
            break;
        end
    else
        if strcmpi(opts.display, 'iter')
            fprintf('Residual did not decrease. Aborting.\n');
        end
        break;
    end

    if kkk == opts.max_iters && strcmpi(opts.display, 'iter')
        fprintf('Reached max iterations. Aborting.\n');
    end
end

solout = sol;
if nargout > 1
    [res, jac] = jacresfunction(solout);
end

% Force solution to be real. TODO: Fix this in cost function instead.
% solout.u = real(solout.u);
% solout.v = real(solout.v);
% solout.a = real(solout.a);
% solout.b = real(solout.b);
% solout.o = real(solout.o);
% 

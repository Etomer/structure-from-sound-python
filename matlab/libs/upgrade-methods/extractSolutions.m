function [Lti,q] = extractSolutions(sample,sols,Hqp,Hqh,opts)
% EXTRACTSOLUTIONS Extracts the L and q from partial upgrade solutions
%   [Lti,q] = EXTRACTSOLUTIONS(sample,sols,Hqp,Hqh,opts) finds the real
%       solutions to L and q provided the solutions from the minimal
%       solvers. Set opts.refine to true if the solutions should be refined
%       using Gauss-Newton. Note: Lti = inv(L').
%
% See also findLinearConstraints, createRandomSample.

    id = sample.id;

    Lti = cell(1,size(sols,2));
    q = zeros(3,size(sols,2));
    goodSol = false(1,size(sols,2));
    for i=1:size(sols,2)
        x = sols(1:id(2)+id(3),i);
        
        if ~all(isfinite(x))
            continue;
        end

        if any(abs(imag(x))./abs(x) > 1e-6)
            continue;
        end
        x = real(x);
        
        if opts.refine
            x = refineUpgrade(sample,Hqp,Hqh,x,opts);
        end
        
        Hq = Hqp+Hqh*x;
        
        H = Hq([1 4 5; 4 2 6; 5 6 3]);
        q(:,i) = Hq([7; 8; 9]);
        
        % This is actually significantly faster than wrapping chol with
        % a try-catch statement.
        if any(eig(H) <= 0)
            continue;
        end

        Lti{i} = chol(H);
        goodSol(i) = 1;
    end
    Lti = Lti(goodSol);
    q = q(:,goodSol);
end


function x = refineUpgrade(sample,Hqp,Hqh,x0,opts)
% REFINEUPGRADE Refines a solution from the action matrix solvers
%   x = REFINEUPGRADE(sample,Hqp,Hqh,x0,opts) takes a single solution x0
%       from one of the action matrix solvers and refines it using
%       Gauss-Newton.
%
% See also extractSolutions.

    maxIters = opts.maxIters;
    tol = opts.tol;
    
    id = sample.id;
    V = sample.V;
    a = sample.a;
    c = sample.c;
    ac = a-c;

    indsI = [1 2 3 1 1 2];
    indsJ = [1 2 3 2 3 3];
    mult = [1 1 1 2 2 2]; % Due to symmetry in H.

    % Some indices and coefficients used when constructing the Jacobian.
    adjHInd = [7 3 2 7 7 6; 7 7 4 3 6 5; 7 5 7 6 2 4; 7 7 4 3 6 5;...
        3 7 1 7 5 7; 6 7 7 5 4 1; 7 5 7 6 2 4; 6 7 7 5 4 1; 2 1 7 4 7 7];
    adjHCoeff = [0 1 1 0 0 -2; 0 0 -1 -1 1 1; 0 -1 0 1 -1 1;...
        0 0 -1 -1 1 1; 1 0 1 0 -2 0; -1 0 0 1 1 -1; 0 -1 0 1 -1 1;...
        -1 0 0 1 1 -1; 1 1 0 -2 0 0];    

    x = x0;
%     fprintf('-----------------------------------------------\n');
    for i=1:maxIters
        Hq = Hqp+Hqh*x;

        H = Hq([1 4 5; 4 2 6; 5 6 3]);
        q = Hq([7; 8; 9]);

        adjH = adj(H); % Still symmetric.
        detH = det(H);

        % Create residuals.
        res = zeros(id(2)+id(3),1);
        if id(2)
            % Note: dot(V,adjH*V)' == diag(V'*adjH*V).
            res(1:id(2),1) = dot(V,adjH*V)'+2*V'*adjH*q-detH*ac';
        end
        if id(3)
            res(id(2)+1:id(2)+id(3),1) = q'*adjH*q-detH*c;
        end
        
%         fprintf('res norm=%e\n',norm(res));
        if norm(res) < tol
            break;
        end

        % Create Jacobian.
        diffAdjH = adjHCoeff.*Hq(adjHInd);
        diffDetH = mult.*adjH(sub2ind([3 3],indsI,indsJ)); % Jacobi's formula.

        vq = repmat(V,3,1).*(repelem(V,3,1)+2*repelem(q,3,1));
        qq = q*q';
        
        jacHb = zeros(id(2)+id(3),9);
        if id(2)
            jacHb(1:id(2),:) = [vq(:,1)'*diffAdjH-ac'*diffDetH 2*V'*adjH];
        end
        if id(3)
            jacHb(id(2)+1:id(2)+id(3),:) = [qq(:)'*diffAdjH-diffDetH*c 2*q'*adjH];
        end
        jacX = jacHb*Hqh;

        % Gauss-Newton step.
        x = x-jacX\res; % TODO: Maybe add some regularization? Levenberg-Marquardt?
    end
end


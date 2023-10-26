function [Lti,q] = upgradeGeneral(sample,opts,actionMatrixSolver)
% UPGRADEGENERAL Solve a low-rank affine upgrade problem
%   [Lti,q] = UPGRADEGENERAL(sample,opts,actionMatrixSolver) finds the
%       affine upgrade L and q given a minimal sample and matching minimal
%       solver. Note that Lti = inv(L').
%
% See also upgrade900, upgrade810, upgrade801.

    id = sample.id;
    V = sample.V;
    a = sample.a;
    c = sample.c;

    [Hqp,Hqh] = findLinearConstraints(sample);
    
    ac = a-c;
    
    ac2 = ac(1:id(2));
    V2 = V(:,1:id(2));
    c2 = c(1:id(3));
    
    sols = actionMatrixSolver([Hqp; Hqh(:); V2(:); ac2(:); c2]);
    
    [Lti,q] = extractSolutions(sample,sols,Hqp,Hqh,opts);
end


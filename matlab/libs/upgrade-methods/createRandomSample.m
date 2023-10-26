function [sample,newProblem] = createRandomSample(id,problem)
% CREATERANDOMSAMPLE Create a random sample from a low-rank upgrade problem
%   [sample,newProblem] = CREATERANDOMSAMPLE(id,problem) creates a random
%       minimal sample of data from a low-rank upgrade problem. The
%       returned sample can be passed to the minimal solver specified by
%       id. The sample is not simply a subset of the full low-rank data
%       (U,V,a,b,c). It is also transformed in a certain way (see
%       Section 3.1. in the paper). newProblem returns the full low-rank
%       data but with these transformations applied.

    U = problem.U;
    V = problem.V;
    a = problem.a;
    b = problem.b;
    c = problem.c;

    m = size(U,2);
    n = size(V,2);
    
    % Random sample.
    indU = randperm(m,id(1)+1);
    indV = randperm(n,id(2)+1);
    
    % Adjust position of the zero columns in U and V (6 dof).
    u0 = U(:,indU(1));
    v0 = V(:,indV(1));
    
    newU = U-u0;
    newV = V-v0;
    newa = a-2*u0'*V;
    newb = b-2*U'*v0;
    newc = c-2*u0'*v0;
    
    % Make sure a(indV(1)) = b(indU(1)) = c (1 dof).
    % This is due to a*ws = wr'*b = c.
    offset = newa(indV(1))-newb(indU(1));
    newa = newa-offset/2;
    newb = newb+offset/2;
    offset = newc-newb(indU(1));
    newa = newa-offset;
    newb = newb-offset;
    newc = newc-offset*2;
    
    sample = struct();
    sample.id = id;
    sample.U = newU(:,indU(2:end));
    sample.V = newV(:,indV(2:end));
    sample.a = newa(indV(2:end));
    sample.b = newb(indU(2:end));
    sample.c = newc;
    sample.indU = indU;
    sample.indV = indV;
    sample.fullU = newU;
    sample.fullV = newV;
    
    newProblem = problem;
    newProblem.U = newU;
    newProblem.V = newV;
    newProblem.a = newa;
    newProblem.b = newb;
    newProblem.c = newc;
end


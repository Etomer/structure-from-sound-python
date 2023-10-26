function [R,t] = kabsch(p1,p2)
%REG_KABSCH Find rigid transformation
%   [R,t] = kabsch(p1,p2) finds the rigid transformation (R,t) that maps p1
%       onto p2, i.e., p2 = R*p1+t.
    
    % Ignore NaNs.
    ind = all(isfinite(p1) & isfinite(p2),1);
    
    p1 = p1(:,ind);
    p2 = p2(:,ind);

    c1 = mean(p1,2);
    c2 = mean(p2,2);

    q1 = p1-c1;
    q2 = p2-c2;

    [U,~,V] = svd(q2*q1');

    R = U*V';
    t = c2-mean(R*p1,2);
end
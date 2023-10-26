function [Lti,q] = upgrade900(sample,~)
% UPGRADE900 Solve a low-rank affine upgrade problem
%   [Lti,q] = UPGRADE900(sample,opts) finds the affine upgrade L and q
%       given a minimal sample with id 900. This solver is linear. Note
%       that Lti = inv(L').
%
% See also upgradeGeneral, upgrade810, upgrade801.

    Hq = findLinearConstraints(sample);
    
    H = Hq([1 4 5; 4 2 6; 5 6 3]);
    q = Hq([7; 8; 9]);

    try
        Lti{1} = chol(H);
    catch
        Lti = {};
        q = [];
    end
end


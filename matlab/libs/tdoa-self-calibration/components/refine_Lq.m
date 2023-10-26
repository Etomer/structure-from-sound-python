function [sol, res, jac] = refine_Lq(sol)
% REFINE_LQ Perform local optimization over upgrade parameters
%   sol = REFINE_LQ(sol) performs local optimization over the upgrade
%       parameters L and q. L and q takes the relaxed solution in u, v, a,
%       b and o to a solution in r, s, and o. Note that this is
%       significantly fewer parameters than when optimizing over the
%       receiver and sender positions.

d = sol.z(sol.rows,sol.cols) - sol.o;
u = sol.u';
v = sol.v;
L = sol.L;
q = sol.q;
inliers = sol.inlmatrix(sol.rows,sol.cols);

[I, J] = find(inliers);
ind = sub2ind(size(d), I, J);
D = d(ind);

[Lopt, qopt, res, jac] = bundletoa_calib(D, I, J, u, v, L, q);

sol.L = Lopt;
sol.q = qopt;
sol.r = Lopt' \ sol.u';
sol.s = Lopt * (sol.v + qopt);

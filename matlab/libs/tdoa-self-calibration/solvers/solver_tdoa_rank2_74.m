function sols = solver_tdoa_rank2_74(d)
% SOLVER_TDOA_RANK2_74 Solve TDOA offsets for the minimal case (7r,4s)
%   sols = SOLVER_TDOA_RANK2_74(z) solves for the TDOA offsets o, given an
%       n=7 by m=4 matrix of TDOA measurements z. The measurements are
%       given by z_ij = || r_i - s_j || + o_j, where r_i and s_j are
%       receiver and senders positions in 2D.

d2 = (d).^2;
u = [(d2(:,2:4)-repmat(d2(:,1),1,3)) (-2*d(:,2:end)) (2*d(:,1))]\ones(7,1);
sols = [u(end)/sum(u(1:3));u((1:3)+3)./u(1:3)];

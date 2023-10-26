function sols = solver_tdoa_rank3_95(z)
% SOLVER_TDOA_RANK3_95 Solve TDOA offsets for the minimal case (9r,5s)
%   sols = SOLVER_TDOA_RANK3_95(z) solves for the TDOA offsets o, given an
%       n=9 by m=5 matrix of TDOA measurements z. The measurements are
%       given by z_ij = || r_i - s_j || + o_j, where r_i and s_j are
%       receiver and senders positions in 3D.

d2 = (z).^2;
u = [(d2(:,2:5)-repmat(d2(:,1),1,4)) (-2*z(:,2:end)) (2*z(:,1))]\ones(9,1);
sols = [u(end)/sum(u(1:4));u((1:4)+4)./u(1:4)];

function sol = init_rso_random(z)
% INIT_RSO_RANDOM Initialize a TDOA solution randomly
%   sol = INIT_RSO_RANDOM(z) initializes a solution in receivers, senders
%       and offsets randomly based on the provided TDOA measurements.

[m, n] = size(z);

sol.rows = 1:m;
sol.cols = 1:n;
sol.inlmatrix = isfinite(z);
zrange = std(z(sol.inlmatrix));
sol.z = z;
sol.r = randn(3, m) * zrange;
sol.s = randn(3, n) * zrange;
sol.o = randn(1, n) * zrange;
sol.type = 'rso';

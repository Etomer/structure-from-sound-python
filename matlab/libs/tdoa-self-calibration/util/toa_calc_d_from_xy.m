function d = toa_calc_d_from_xy(x,y)
% TOA_CALC_D_FROM_XY Calculate pair-wise distances between points
% d = TOA_CALC_D_FROM_XY(x,y) calculates the pairwise Euclidean distance
%   between the points in x and y.
%   Input: 
%       d - m x n matrix
%   Output: 
%       x - d x m matrix
%       y - d x n matrix

[x_dim,m] = size(x);
[y_dim,n] = size(y);

if x_dim > y_dim
    y((y_dim+1):x_dim,:) = zeros(x_dim-y_dim,n);
elseif y_dim > x_dim
    x((x_dim+1):y_dim,:) = zeros(y_dim-x_dim,n);
end

d = sqrt(sum((kron(ones(1,n),x) - kron(y,ones(1,m))).^2 , 1));
d = reshape(d,m,n);

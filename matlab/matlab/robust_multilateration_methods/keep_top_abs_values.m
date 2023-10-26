function M_out = keep_top_abs_values(M_in,K);
% M_out = keep_top_abs_values(M_in,K);

[m,n]=size(M_in);
mv = M_in(:);
[sortv,sorti]=sort(abs(mv));
mv(sorti(1:(end-K))) = zeros(length(mv)-K,1);
M_out = reshape(mv,m,n);

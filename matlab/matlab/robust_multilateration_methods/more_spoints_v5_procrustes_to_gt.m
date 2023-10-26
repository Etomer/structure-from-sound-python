function asol_out = more_spoints_v5_procrustes_to_gt(asol,rgt);


r = asol.r;
s = asol.s;

mean(sqrt(sum( (rgt-r).^2 )))
[D,Z,T] = procrustes(rgt',r')
r_new = T.b * r' * T.T + repmat(T.c(1,:),size(r,2),1);
r_new = r_new';
s_new = T.b * s' * T.T + repmat(T.c(1,:),size(s,2),1);
s_new = s_new';
mean(sqrt(sum( (rgt-r_new).^2 )))

asol_out = asol;
asol_out.r = r_new;
asol_out.s = s_new;


function [res,jac]=misstdoa_jacres_v2_fixneg(sol);


C = 1;

if 0,
    sol0 = sol;
    kk = 1; dz = zeros(128,1); dz(kk)=1;
    litet = 0.00001;
    sol1 = misstdoa_update_v2(sol,litet*dz);
    [res0,jac0]=misstdoa_jacres_v2(sol0);
    [res1,jac1]=misstdoa_jacres_v2(sol1);
    [(res1-res0)/litet jac0*dz]
  
    sol0 = sol;
    dz = randn(128,1); dz = dz/norm(dz);
    litet = 0.00001;
    sol1 = misstdoa_update_v2(sol,litet*dz);
    [res0,jac0]=misstdoa_jacres_v2(sol0);
    [res1,jac1]=misstdoa_jacres_v2(sol1);
    [(res1-res0)/litet jac0*dz]

    for kk = 1:128,
        sol0 = sol;
        dz = zeros(128,1); dz(kk)=1;
        litet = 0.00001;
        sol1 = misstdoa_update_v2(sol,litet*dz);
        [res0,jac0]=misstdoa_jacres_v2(sol0);
        [res1,jac1]=misstdoa_jacres_v2(sol1);
        [(res1-res0)/litet jac0*dz]
        kk
        pause;
        
    end;
    
    for kk = 1:128,
        sol0 = sol;
        dz = zeros(128,1); dz(kk)=1;
        litet = 0.00001;
        sol1 = misstdoa_update_v2(sol,litet*dz);
        [res0,jac0]=misstdoa_jacres_v2(sol0);
        [res1,jac1]=misstdoa_jacres_v2(sol1);
        [(res1-res0)/litet jac0*dz]
        [norm((res1-res0)/litet) norm((res1-res0)/litet-jac0*dz)]
        kk
        
        pause;
        
    end;
    
    
end;

z = sol.z;
zcut = z(sol.rows,sol.cols);
icut = sol.inlmatrix(sol.rows,sol.cols);
[I,J]=find(icut);
zind = sub2ind(size(zcut),I,J);
z_vect = zcut(zind); % vector of all z measurements that are inliers

m1 = length(sol.rows);
n1 = length(sol.cols);
M = length(I);
% residual is size Mx1
% parameter order
% w1 w2 w3 v1 v2 v3 c d o
N = 3*m1+3*n1+m1+n1+n1;
% jacobean is size MxN
%
% internal
in = -2*sol.u(I,1).*(sol.v(1,J)')  -2*sol.u(I,2).*(sol.v(2,J)')-2*sol.u(I,3).*(sol.v(3,J)')+ sol.a(I)+sol.b(J)';
% residuals
%res0 = sqrt(in)+sol.o(J)'-z_vect;
res1 = sqrt(relu(in))+sol.o(J)'-z_vect;
res2 = relu(-in);

% jacobean
sin = sqrt(relu(in));
%ds = (1/2)*1./sin;
ds = zeros(size(in));
posid = find(in>0);
ds(posid)=(1/2)*1./sin(posid);
Jw1 = 1:m1;
Jw2 = (1:m1)+m1;
Jw3 = (1:m1)+2*m1;
Jv1 = (1:n1) + 3*m1;
Jv2 = (1:n1)+3*m1+1*n1;
Jv3 = (1:n1)+3*m1+2*n1;
Jc  = (1:m1)+3*m1+3*n1;
Jd  = (1:n1)+4*m1+3*n1;
Jo  = (1:n1)+3*m1+4*n1;
%
II = repmat((1:M)',9,1);
Jw1 = I;
Jw2 = I+m1;
Jw3 = I+2*m1;
Jv1 = J + 3*m1;
Jv2 = J+3*m1+1*n1;
Jv3 = J+3*m1+2*n1;
Jc  = I+3*m1+3*n1;
Jd  = J+4*m1+3*n1;
Jo  = J+4*m1+4*n1;
JJ = [Jw1;Jw2;Jw3;Jv1;Jv2;Jv3;Jc;Jd;Jo];
%
dw1 = -2*ds.*(sol.v(1,J)');
dw2 = -2*ds.*(sol.v(2,J)');
dw3 = -2*ds.*(sol.v(3,J)');
dv1 = -2*ds.*(sol.u(I,1));
dv2 = -2*ds.*(sol.u(I,2));
dv3 = -2*ds.*(sol.u(I,3));
dc  = ds;
dd  = ds;
do  = ones(M,1);
DD = [dw1;dw2;dw3;dv1;dv2;dv3;dc;dd;do];
%jac1 = sparse(II,JJ,DD,M,N);


%in = -2*sol.u(I,1).*(sol.v(1,J)')  -2*sol.u(I,2).*(sol.v(2,J)')-2*sol.u(I,3).*(sol.v(3,J)')+ sol.a(I)+sol.b(J)';
II2 = repmat((1:M)'+M,9,1);
dw12 = -2*(sol.v(1,J)');
dw22 = -2*(sol.v(2,J)');
dw32 = -2*(sol.v(3,J)');
dv12 = -2*(sol.u(I,1));
dv22 = -2*(sol.u(I,2));
dv32 = -2*(sol.u(I,3));
dc2  = ones(M,1);
dd2  = ones(M,1);
do2 = zeros(M,1);
DD2 = C*[dw12;dw22;dw32;dv12;dv22;dv32;dc2;dd2;do2];
%jac2 = sparse(II2,JJ,DD2,2*M,N);

%%
res = [res1;C*res2];
jac = sparse([II;II2],[JJ;JJ],[DD;DD2],2*M,N);
% res = [res1];
% jac = sparse([II],[JJ],[DD],M,N);



function [y,res,jac]=tdoa_multilaterate_truncl1(data,x,y0);
% y=toa_trilateration_one_bundle(d,x,y0)
% non-linear least squares optimization of y with y0
% as initial estimate, i e
% minimise
%  min_y  sum_i (d_i - sqrt(sum( (x(:,i)-y).^2 )))

if 0,
    xt = x;
    yt = y0;
    dz = [1;0;0];
    %   dz = [0;1;0];
    dz = [0;0;1];
    litet = 0.0001;
    y0 = yt;
    [y1]=updatey(y0,dz*litet);
    [res0,jac0]=calcresandjac(data,xt,y0);
    [res1,jac1]=calcresandjac(data,xt,y1);
    [(res1-res0)/litet jac0*dz]
end





debug = 0;
%    keyboard;
% Initial weights
w = ones(1,size(data,2))';
e = w;
% Minimize sum( w.*(res.^2) );
xt = x;
yt = y0;
%keyboard;


for kkk = 1:20;
    %kkk
    [res,jac]=calcresandjac(data,xt,yt);
    w = (abs(res) <0.2).* sqrt((e ./ (abs(res)+10^(-10)))) + ...
        sqrt(0.2)*(abs(res)>=0.2).*     ((e ./ (abs(res)+10^(-10))));
    mres = (w.*res);
    mjac = (repmat(w,1,size(jac,2)).*jac);
    dz = -(mjac\mres);
    %dz = -(jac'*jac+eye(size(jac,2)))\(jac'*res);
    %[u,s,v]=svd(full(jac),0);
    %keyboard;
    %     nrparam = size(jac,2);
    %     dof = nrparam-6;
    %     u = u(:,1:dof);
    %     s = s(1:dof,1:dof);
    %     v = v(:,1:dof);
    [ytn]=updatey(yt,dz);
    [res2,jac2]=calcresandjac(data,xt,ytn);
    %w2 = sqrt((e ./ (abs(res)+10^(-10))));
    mres2 = (w.*res2);
    aa = [norm(res) norm(res+jac*dz) norm(res2)];
    aa = [sum(abs(res)) sum(abs(res+jac*dz)) sum(abs(res2)) ];
    aa = [norm(mres) norm(mres+mjac*dz) norm(mres2)];
    bb = aa;
    bb=bb-bb(2);
    bb = bb/bb(1);
    cc = norm(jac*dz)/norm(res);
    % Check that the error actually gets smaller
    if norm(mres)<norm(mres2),
        % bad
        % check that the update is sufficiently big
        % otherwise it is maybe just numerical limitations
        if cc>1e-4,
            % OK then it is probably just the linearization that
            % is not good enough for this large step size, decrease
            kkkk = 1;
            while (kkkk<50) & (norm(mres)<norm(mres2)),
                dz = dz/2;
                [ytn]=updatey(yt,dz);
                [res2,jac2]=calcresandjac(data,xt,ytn);
                mres2 = (w.*res2);
                kkkk = kkkk+1;
            end
        end
    end
    if debug,
        aa = [norm(mres) norm(mres+mjac*dz) norm(mres2)];
        bb = aa;
        bb=bb-bb(2);
        bb = bb/bb(1);
        cc = norm(jac*dz)/norm(res);
        %keyboard;
        aa
        bb
        cc
    end;
    if norm(mres2)<norm(mres)
        yt = ytn;
        res = res2;
    else
        %disp([num2str(kkk) '  stalled']);
    end
    %sum(abs(res))
end;

y = yt;


function [res,jac]=calcresandjac(data,x,y0);

%keyboard;
nn = size(data,2);
I1 = data(1,:);
I2 = data(2,:);
d = data(3,:)';
V1 = x(:,I1)-repmat(y0,1,nn);
V2 = x(:,I2)-repmat(y0,1,nn);
Vt1 = V1';
Vt2 = V2';
dd1 = sqrt(sum(V1.^2,1))';
dd2 = sqrt(sum(V2.^2,1))';
res = dd2-dd1-d;
idd1 = 1./dd1;
idd2 = 1./dd2;
II = (1:nn)';
JJ1 = 1*ones(size(II));
JJ2 = 2*ones(size(II));
JJ3 = 3*ones(size(II));
VV1 = -idd2.*Vt2(:,1)+idd1.*Vt1(:,1);
VV2 = -idd2.*Vt2(:,2)+idd1.*Vt1(:,2);
VV3 = -idd2.*Vt2(:,3)+idd1.*Vt1(:,3);
jac = sparse([II;II;II],[JJ1;JJ2;JJ3],[VV1;VV2;VV3],nn,3);


function [yny]=updatey(y,dz);
yny = y + dz(1:end);




function [Lopt,qopt,res,jac]=bundletoa_calib(D,I,J,u,v,Lt,qt,debug);
% [xopt,yopt,res,jac]=bundletoa_calib(D,I,J,ut,vt,Lt,qt,debug);
%

if nargin<8,
    debug = 0;
end;

if 0,
    % Tests to check derivatives
    [L0,q0]=updateLq(Lt,qt,0.1*randn(9,1));
    for j = 1:9,
        dz = zeros(9,1); dz(j) = 1;
        litet = 0.00001;
        [L1,q1]=updateLq(L0,q0,litet*dz);
        [res0,jac0]=calcresandjac2(D,I,J,u,v,L0,q0);
        [res1,jac1]=calcresandjac2(D,I,J,u,v,L1,q1);
        [(res1-res0)/litet jac0*dz]
        [j norm((res1-res0)/litet) norm(jac0*dz) norm( (res1-res0)/litet-jac0*dz )]
        pause;
    end
    for j = 1:9,
        dz = randn(9,1); dz = dz/norm(dz);
        litet = 0.00001;
        [L1,q1]=updateLq(L0,q0,litet*dz);
        [res0,jac0]=calcresandjac2(D,I,J,u,v,L0,q0);
        [res1,jac1]=calcresandjac2(D,I,J,u,v,L1,q1);
        [(res1-res0)/litet jac0*dz]
        [j norm((res1-res0)/litet) norm(jac0*dz) norm( (res1-res0)/litet-jac0*dz )]
        pause;
    end
end


for kkk = 1:30;
    %kkk
    %keyboard;
    [res,jac]=calcresandjac2(D,I,J,u,v,Lt,qt);
    %dz = -(jac\res);
    %dz = -(jac'*jac+speye(size(jac,2)))\(jac'*res);
    dz = -(jac'*jac+0.001*speye(size(jac,2)))\(jac'*res);
    %[u,s,v]=svd(full(jac),0);
    %keyboard;
    %nrparam = size(jac,2);
    %dof = nrparam-6;
    %u = u(:,1:dof);
    %s = s(1:dof,1:dof);
    %v = v(:,1:dof);
    %dz = -v*inv(s)*u'*res;
    [Ltn,qtn]=updateLq(Lt,qt,dz);
    [res2,jac2]=calcresandjac2(D,I,J,u,v,Ltn,qtn);
    aa = [norm(res) norm(res+jac*dz) norm(res2)];
    bb = aa;
    bb=bb-bb(2);
    bb = bb/bb(1);
    cc = norm(jac*dz)/norm(res);
    % Check that the error actually gets smaller
    if norm(res)<norm(res2),
        % bad
        % check that the update is sufficiently big
        % otherwise it is maybe just numerical limitations
        if cc>1e-4,
            % OK then it is probably just the linearization that
            % is not good enough for this large step size, decrease
            kkkk = 1;
            while (kkkk<50) & (norm(res)<norm(res2)),
                dz = dz/2;
                [Ltn,qtn]=updateLq(Lt,qt,dz);
                [res2,jac2]=calcresandjac2(D,I,J,u,v,Ltn,qtn);
                kkkk = kkkk+1;
            end
        end
    end
    if debug,
        aa = [norm(res) norm(res+jac*dz) norm(res2)];
        bb = aa;
        bb=bb-bb(2);
        bb = bb/bb(1);
        cc = norm(jac*dz)/norm(res);
        %keyboard;
        aa
        bb
        cc
    end;
    if norm(res2)<norm(res)
        Lt = Ltn;
        qt = qtn;
    else
        %disp([num2str(kkk) '  stalled']);
    end
end;

Lopt = Lt;
qopt = qt;

function [res,jac]=calcresandjac2(D,I,J,u,v,L,q);

m = size(u,2);
n = size(v,2);
% calculate upgraded positions (x,y) from (u,v) using L and q
x = L*(u+repmat(q,1,size(u,2)));
y = inv(L')*v;
% calculated res and jac with respect to changes in x,y
[res,jac_xy]=calcresandjac(D,I,J,x,y);
% calculate d(x,y)/d(L,q)
Lbase{1} = [1 0 0;0 0 0;0 0 0];
Lbase{2} = [0 1 0;0 0 0;0 0 0];
Lbase{3} = [0 0 1;0 0 0;0 0 0];
Lbase{4} = [0 0 0;0 1 0;0 0 0];
Lbase{5} = [0 0 0;0 0 1;0 0 0];
Lbase{6} = [0 0 0;0 0 0;0 0 1];
qbase{1} = [1 0 0]';
qbase{2} = [0 1 0]';
qbase{3} = [0 0 1]';
dxydLq = zeros(3*m+3*n,9);
for j = 1:6
    dxdL = Lbase{j}*(u+repmat(q,1,size(u,2)));
    dydL = -inv(L')*Lbase{j}'*inv(L')*v;
    dxydLq(:,j) = [dxdL(:);dydL(:)];
end
for j = 1:3
    dxdb = L*(repmat(qbase{j},1,size(u,2)));
    dydb = zeros(3,n);
    dxydLq(:,6+j) = [dxdb(:);dydb(:)];
end
% Use chain-rule to get jac with respect to changes in L,q
jac = jac_xy*dxydLq;


function [res,jac]=calcresandjac(D,I,J,x,y);

nn = length(D);
m = size(x,2);
n = size(y,2);
V = x(:,I)-y(:,J);
Vt = V';
dd = sqrt(sum(V.^2,1))';
idd = 1./dd;
res = dd-D;
II = (1:length(I))';
JJ1 = (I-1)*3+1;
JJ2 = (I-1)*3+2;
JJ3 = (I-1)*3+3;
JJ4 = (J-1)*3+1+3*m;
JJ5 = (J-1)*3+2+3*m;
JJ6 = (J-1)*3+3+3*m;

VV1 = idd.*Vt(:,1);
VV2 = idd.*Vt(:,2);
VV3 = idd.*Vt(:,3);
VV4 = -idd.*Vt(:,1);
VV5 = -idd.*Vt(:,2);
VV6 = -idd.*Vt(:,3);

jac = sparse([II;II;II;II;II;II],[JJ1;JJ2;JJ3;JJ4;JJ5;JJ6],[VV1;VV2;VV3;VV4;VV5;VV6],nn,3*m+3*n);



function [Lny,qny]=updateLq(L,q,dz);

Lny = L;
qny = q;
Lbase{1} = [1 0 0;0 0 0;0 0 0];
Lbase{2} = [0 1 0;0 0 0;0 0 0];
Lbase{3} = [0 0 1;0 0 0;0 0 0];
Lbase{4} = [0 0 0;0 1 0;0 0 0];
Lbase{5} = [0 0 0;0 0 1;0 0 0];
Lbase{6} = [0 0 0;0 0 0;0 0 1];
qbase{1} = [1 0 0]';
qbase{2} = [0 1 0]';
qbase{3} = [0 0 1]';
for j = 1:6
    Lny = Lny+Lbase{j}*dz(j);
end
for j = 1:3
    qny = qny+qbase{j}*dz(6+j);
end





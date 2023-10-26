function asol = more_spoints_v0_just_put_in_firssol(detections,firstsol,u,exp);

% First make a new representation

tn = size(detections.toas,2);
stmp = NaN*ones(3,tn);
stmp(:,firstsol.okcols)=firstsol.s;
% % Hur ska man representera mätningarna
% i1 = mikrofon 1
% i2 = mikrofon 2
% j = tidpunkt
% k = vilken topp av de fyra bästa från gcc-phat?
% om vi vill kunna ha speglade mikrofoner så måste vi 
% ha en mappning till vilken ljudkanal vi ska leta efter
m2l = 1:12;
I1 = zeros(1,tn*12);
I2 = zeros(1,tn*12);
J = zeros(1,tn*12);
K = zeros(1,tn*12);
U = zeros(1,tn*12);
UC = zeros(1,tn*12);
kk = 0;
%
for jk = 1:length(firstsol.okcols),
    %jk = 3134;
    j = firstsol.okcols(jk);
    ss = firstsol.s(:,jk);
    i1 = 6;
    for i2 = setdiff(1:12,i1),
        uc_tmp = detections.z(i2,j);
        [minv,mink]=min(abs(u{i1,i2}(:,j)/96000*exp.speedofsound-uc_tmp));
        k = mink;
        u_tmp = norm(firstsol.r(:,i2)-ss)-norm(firstsol.r(:,i1)-ss);
        if abs(u_tmp-uc_tmp)<0.1,
            kk = kk+1;
            I1(kk)=i1;
            I2(kk)=i2;
            J(kk)=j;
            K(kk)=k;
            U(kk)=u_tmp;
            UC(kk)=uc_tmp;
        end
    end;
end
I1 = I1(1:kk);
I2 = I2(1:kk);
J = J(1:kk);
K = K(1:kk);
U = U(1:kk);
UC = UC(1:kk);

asol.r = firstsol.r;
asol.s = stmp;
asol.I1 = I1;
asol.I2 = I2;
asol.J = J;
asol.K = K;
asol.U = U;
asol.UC = UC;


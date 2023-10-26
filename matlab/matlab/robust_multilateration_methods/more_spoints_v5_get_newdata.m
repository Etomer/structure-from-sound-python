function asol_out = more_spoints_v5_get_newdata(asol,u,speedofsound,a_sr,threshold_dist);

if nargin<5,
    threshold_dist = 0.1;
end;

threshold_dist

threshold_nr = 0;
s_in = asol.s;
r = asol.r;
nrinls_out = zeros(1,size(s_in,2));

jj = 1:size(s_in,2);
jj = jj(find(isfinite(s_in(1,:))));

I1 = [];
I2 = [];
J = [];
K = [];
U = [];
UC = [];
nrinls_out = zeros(1,size(s_in,2));
%keyboard;
%
for j = 1:size(s_in,2),
    
    if round(j/100)==(j/100)
        disp([j size(s_in,2)]);
    end

    
    if 1, % Get the data from nearby time instances
        jjj = j;
        udata = [];
        for i1 = 1:11,
            for i2 = (i1+1):12,
                udata = [udata [i1*ones(1,4);i2*ones(1,4);u{i1,i2}(:,jjj)';(jjj-j)*ones(1,4)]];
            end
        end
    end;
    behall = find(isfinite(udata(3,:)));
    udata = udata(:,behall);
    behall = find(udata(3,:)~=0);
    udata = udata(:,behall);
    udata(3,:)=udata(3,:)*speedofsound/a_sr;
    
    %jk = 3134;
    etts = s_in(:,j);
    
    if isfinite(etts(1,1)),
        reproj = sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(2,:))).^2 )) - ...
            sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(1,:))).^2 ));
        res = reproj-udata(3,:);
        nrinl = sum( abs(res)<threshold_dist );
        if nrinl>= threshold_nr,
            inl = find(abs(res)<threshold_dist);
            udata = udata(:,inl);
            reproj = reproj(:,inl);
            nrinls_out(1,j)=nrinl;
            I1 = [I1 udata(1,:)];
            I2 = [I2 udata(2,:)];
            J  = [J j*ones(1,size(udata,2))];
            K = [K 1*ones(1,size(udata,2))];
            U = [U udata(3,:)];
            UC = [UC reproj];
        end
    end
end;


asol_out.nrinls = nrinls_out;
asol_out.r = asol.r;
asol_out.s = asol.s;
asol_out.I1 = I1;
asol_out.I2 = I2;
asol_out.J = J;
asol_out.K = K;
asol_out.U = U;
asol_out.UC = UC;

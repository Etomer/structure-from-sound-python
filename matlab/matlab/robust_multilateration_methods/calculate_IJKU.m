function asol_out = calculate_IJKU(asol,u,speedofsound,a_sr);

threshold = 0.1;
inl_threshold = 10;

s_in = asol.s;
s_out = s_in;
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

for j = 1:size(s_out,2);
    
    
    if round(j/1000)==(j/1000)
        fprintf('%s',['-' num2str(j/1000)]);
    end
    
    
    if 1,
        udata = [];
        for i1 = 1:11,
            for i2 = (i1+1):12,
                udata = [udata [i1*ones(1,4);i2*ones(1,4);u{i1,i2}(:,j)']];
            end
        end
    end;
    udata(3,:)=udata(3,:)*speedofsound/a_sr;
    behall = find(isfinite(udata(3,:)));
    udata = udata(:,behall);
    
    
    
    
    etts = s_out(:,j);
    if isfinite(etts(1,1)),
        reproj = sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(2,:))).^2 )) - ...
            sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(1,:))).^2 ));
        res = reproj-udata(3,:);
        nrinl = sum( abs(res)<threshold );
        if nrinl>= inl_threshold,
            inl = find(abs(res)<threshold);
            udata = udata(:,inl);
            reproj = reproj(inl);
            stmp = etts; % Don't do any additional optimization. Just calculate the statistics
            reproj = sqrt(sum( (repmat(stmp,1,size(udata,2))-r(:,udata(2,:))).^2 )) - ...
                sqrt(sum( (repmat(stmp,1,size(udata,2))-r(:,udata(1,:))).^2 ));
            s_out(:,j)=stmp;
            nrinls_out(1,j)=nrinl;
            I1 = [I1 udata(1,:)];
            I2 = [I2 udata(2,:)];
            J  = [J j*ones(1,size(udata,2))];
            K = [K 1*ones(1,size(udata,2))];
            U = [U udata(3,:)];
            UC = [UC reproj];
        else
        end
    end
end;

asol_out.r = asol.r;
asol_out.s = s_out;
asol_out.I1 = I1;
asol_out.I2 = I2;
asol_out.J = J;
asol_out.K = K;
asol_out.U = U;
asol_out.UC = UC;

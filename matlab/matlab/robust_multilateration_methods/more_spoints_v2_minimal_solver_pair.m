function asol_out = more_spoints_v2_minimal_solver_pair(asol,u,speedofsound,a_sr);

threshold = 0.1;
inl_threshold = 10;

s_in = asol.s;
r = asol.r;
s_out = NaN*ones(size(s_in));
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
    
    if size(udata,2)>10,
        % Ta startgissning från minimallösaren
        maxinl = 0;
        for kkk = 1:1000;
            sel = randperm(size(udata,2),3);
            tmp1 = udata(:,sel);
            tmp2 = tmp1(1:2,:);
            tmp3 = r(:,tmp2(:));
            tmp4 = [tmp3(:); tmp1(3,:)'.^2];
            try
                sols = solver_tdoa_pair_3d(tmp4);
                for k = 1:size(sols,2);
                    etts = real(sols(:,k));
                    reproj = sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(2,:))).^2 )) - ...
                        sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(1,:))).^2 ));
                    res = reproj-udata(3,:);
                    [sortv,sorti]=sort(abs(res));
                    if kkk<3,
                        %figure(30); plot(log10(sortv+0.001),'*'); hold on;
                    end
                    nrinl = sum( abs(res)<threshold );
                    
                    if nrinl>maxinl,
                        maxinl = nrinl;
                        disp([j maxinl]);
                        bests = etts;
                        %figure(30); plot(log10(sortv+0.001),'*'); hold on;
                        nrinl = sum( abs(res)<threshold );
                        
                    end
                end
            catch
            end
        end
        
        
        
        etts = bests;
        if isfinite(etts(1,1)),
            reproj = sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(2,:))).^2 )) - ...
                sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(1,:))).^2 ));
            res = reproj-udata(3,:);
            nrinl = sum( abs(res)<threshold );
            if nrinl>= inl_threshold,
                inl = find(abs(res)<threshold);
                udata = udata(:,inl);
                reproj = reproj(inl);
                [stmp,res_out,jac_out]=tdoa_multilaterate(udata,r,etts);
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
end

asol_out.r = asol.r;
asol_out.s = s_out;
asol_out.I1 = I1;
asol_out.I2 = I2;
asol_out.J = J;
asol_out.K = K;
asol_out.U = U;
asol_out.UC = UC;

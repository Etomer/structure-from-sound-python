function asol_out = more_spoints_v3_use_interval_1(asol,u,speedofsound,a_sr);

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
    if round(j/100)==(j/100)
        disp([j size(s_out,2)]);
    end
    if 1, % Get the data from nearby time instances
        udata = [];
        for jjj = max((j-3),1):min((j+3),size(asol.s,2)),
            for i1 = 1:11,
                for i2 = (i1+1):12,
                    udata = [udata [i1*ones(1,4);i2*ones(1,4);u{i1,i2}(:,jjj)';(jjj-j)*ones(1,4)]];
                end
            end
        end;
    end;
    behall = find(isfinite(udata(3,:)));
    udata = udata(:,behall);
    behall = find(udata(3,:)~=0);
    udata = udata(:,behall);
    udata(3,:)=udata(3,:)*speedofsound/a_sr;
    
    if size(udata,2)>10,
        % Ta startgissning från närliggande punkter från asol
        s_guess = s_in(:,max((j-10),1):min((j+10),size(asol.s,2)));
        s_behall = find( sqrt(sum( s_guess.^2 )) <100 );
        s_guess = s_guess(:,s_behall);
        
        maxinl = 0;
        bests = NaN*ones(3,1);
        for kkk = 1:size(s_guess,2);
            
            etts = s_guess(:,kkk);
            %etts = randn(3,1);
            reproj = sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(2,:))).^2 )) - ...
                sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(1,:))).^2 ));
            res = reproj-udata(3,:);
            
            
            %             restmp = sort(abs(res));
            %             figure(7); hold on;
            %             plot(restmp,'r-');
            %             axis([0 1800 0 0.4]);
            
            %             figure(1); clf;
            %             subplot(2,1,1);
            %             hold off;
            %             plot(udata(3,:),'g*');
            %             hold on;
            %             plot(reproj,'bo');
            %             figure(1); subplot(2,1,2);
            %             plot(res,'*');
            
            [sortv,sorti]=sort(abs(res));
            %             if kkk<3,
            %                 %figure(30); plot(log10(sortv+0.001),'*'); hold on;
            %             end
            nrinl = sum( abs(res)<threshold );
            %disp([kkk nrinl]);
            if nrinl>maxinl,
                maxinl = nrinl;
                %disp([j maxinl]);
                bests = etts;
                
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
                %stmp = etts;
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

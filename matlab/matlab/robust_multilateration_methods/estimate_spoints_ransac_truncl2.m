function asol_out = estimate_spoints_ransac_truncl2(r,u,speedofsound,a_sr);

warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:rankDeficientMatrix')

% We really only need the microphone positions
% r
% from asol
% and the measurement
% u
%
% We might not need to store everything I1,I2 and so forth either

s_out = NaN*ones(3,size(u{1,2},2));
threshold = 0.1;
inl_threshold = 10;

for j = 1:size(s_out,2);
    if round(j/100)==(j/100)
        fprintf('%s',['-' num2str(j/100)]);
    end
    
    % Extract data from one time instant -> udata
    udata = [];
    for i1 = 1:11,
        for i2 = (i1+1):12,
            udata = [udata [i1*ones(1,4);i2*ones(1,4);u{i1,i2}(:,j)']];
        end
    end
    udata(3,:)=udata(3,:)*speedofsound/a_sr;
    behall = find(isfinite(udata(3,:)));
    udata = udata(:,behall);
    
    %
    %keyboard;
    
    if size(udata,2)>10,
        %j
        % Ta startgissning från minimallösaren
        bests = NaN*ones(3,1);
        maxinl = 0;
        for kkk = 1:100;
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
                        %disp([j maxinl]);
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
             else
            end
        end
    end;
    
end

asol_out.r = r;
asol_out.s = s_out;

fprintf('%s\n',['-end']);

warning('on','MATLAB:singularMatrix')
warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:rankDeficientMatrix')

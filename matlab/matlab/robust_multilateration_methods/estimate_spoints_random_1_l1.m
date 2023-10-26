function asol_out = estimate_spoints_random_1_l1(r,u,speedofsound,a_sr);

% We really only need the microphone positions
% r 
% from asol
% and the measurement
% u
% 
% We might not need to store everything I1,I2 and so forth either

s_out = NaN*ones(3,size(u{1,2},2));

for j = 1:size(s_out,2);
    if round(j/1000)==(j/1000)
        fprintf('%s',['-' num2str(j/1000)]);
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
    
    % Use one state-of-the-art method to estimate s
    
    inl_threshold = 3;
    nrinl = size(udata,2);
    if nrinl>= inl_threshold, % We need at least 3 measurements
        
        etts = 5*randn(3,1); % Random initial guess
        [stmp,res_out,jac_out]=tdoa_multilaterate_l1(udata,r,etts);
%         reproj = sqrt(sum( (repmat(stmp,1,size(udata,2))-r(:,udata(2,:))).^2 )) - ...
%             sqrt(sum( (repmat(stmp,1,size(udata,2))-r(:,udata(1,:))).^2 ));
        s_out(:,j)=stmp;
    end
end

asol_out.r = r;
asol_out.s = s_out;

fprintf('%s\n',['-end']);

function asol_out = more_spoints_nonlinearleastsquares(asol,u,speedofsound,a_sr);

% We really only need the microphone positions
% r 
% from asol
% and the measurement
% u
% 
% We might not need to store everything I1,I2 and so forth either

r = asol.r;
s_out = NaN*ones(size(asol.s));

% Behövs nrinls_out och I1,I2,J,K,U,UC? 
% Jag hade det i andra skript, men det är kanske inte nödvändigt. 
nrinls_out = zeros(1,size(s_out,2));  
I1 = [];
I2 = [];
J = [];
K = [];
U = [];
UC = [];

for j = 1:size(s_out,2);
    if round(j/100)==(j/100)
        display([j size(s_out,2)]);
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
    end
end

asol_out.r = asol.r;
asol_out.s = s_out;
asol_out.I1 = I1;
asol_out.I2 = I2;
asol_out.J = J;
asol_out.K = K;
asol_out.U = U;
asol_out.UC = UC;



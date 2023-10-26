function asol_out = estimate_spoints_velasco(r,u,speedofsound,a_sr);

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
    % This time take only the top peak, since Velasco
    % can only handle one peak
    % Save in matrix M instead
    M = zeros(size(u)); % Hardcoded
    for i1 = 1:size(u,1),
        for i2 = 1:size(u,2),
            M(i1,i2)=u{i1,i2}(1,j);
        end
    end
    L = double(isfinite(M));
    for i1 = 1:size(u,1);
        M(i1,i1)=0;
    end
    for i1 = 1:size(u,1);
        for i2 = 1:size(u,2);
            if ~isfinite(M(i1,i2)),
                M(i1,i2)=0;
            end
        end
    end
    
    % Use one state-of-the-art method to estimate s
    
    nrinl = (sum(L(:))-12)/2;
    inl_threshold = 20;
    Db = eye(size(L)).*(L*L');
    E = ones(size(M));
    Lbar = E-L;
    if (nrinl>= inl_threshold) & (abs(det(Db+Lbar))>0.0001), % We need at least 20 measurements
        
        % Use Velasco. Go from M to v
        % algorithm
        Mtilde = (M-M')/2;
        Mtilde = Mtilde*speedofsound/a_sr;
        Mtilde = Mtilde.*L;
        E = ones(size(Mtilde));
        k = 5;
        epsilon = 10^(-10);
        % alg
        Db = eye(size(L)).*(L*L');
        Lbar = E-L;
        Q = inv(Db+Lbar);
        M0 =zeros(size(Mtilde));
        Mt = M0;
        S0 =zeros(size(Mtilde));
        St = S0;
        t=0;
        while (norm(L.*(Mtilde - Mt)-St,'fro')^2/norm(Mtilde,'fro')^2 > epsilon) & (t<20),
            t=t+1;
            Mt = (Q*(L.*Mtilde-St)*E+E*(L.*Mtilde-St)*Q);
            St = keep_top_abs_values(L.*(Mtilde-Mt),2*k);
            norm(L.*(Mtilde - Mt)-St,'fro')^2/norm(Mtilde,'fro')^2;
        end;
        v = -Mt(:,1);
        
        if all(isfinite(v)),
            % Use v to estimate s from r and v
            
            %[ysols,osols]=tdoa_trilateration_y_one_point(v,r);
            %[ysols,osols]=tdoa_trilateration_y_one_point_Chen_Ho(v,r);
            nrok = sum(isfinite(v));
            if nrok >=4,
                [ysol,osol,inliers_ut,nr_inliers_ut,err_rms]=tdoa_trilateration_y_one_ransac(v,r,30,0.10);
                if nr_inliers_ut>=4,
                    [ysol2,osol2,res2]=tdoa_trilateration_y_one_bundle(v,r,ysol,osol,inliers_ut);
                    s_out(:,j)=ysol2;
                end
            end;
        end
    end
end

asol_out.r = r;
asol_out.s = s_out;


fprintf('%s\n',['-end']);
function asol_out = estimate_spoints_ChanHo(r,u,speedofsound,a_sr);

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
    L = isfinite(M);
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
    inl_threshold = 10;
    if nrinl>= inl_threshold, % We need at least 20 measurements
        
        % Use simple method to go from M to v
        v = -M(:,6)*speedofsound/a_sr;
        
        % Use v to estimate s from r and v
        
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
        
        s_out(:,j)=ysol2;
    end
end

asol_out.r = r;
asol_out.s = s_out;


fprintf('%s\n',['-end']);


function [zcalc,zok] = misstdoa_calc_z(sol)
% MISSTDOA_CALC_Z

z = sol.z;
m = length(sol.rows);
n = length(sol.cols);

switch lower(sol.type)
    case {'uvabo'}
        tmp = -2*(sol.u*sol.v) + repmat(sol.a,1,n) + repmat(sol.b,m,1);
        % This is (kind of) -2* U'*V + a_i + b_j
        ok = tmp >= 0;
        neg = tmp < 0;
        ztmp = sqrt(relu(tmp)) + repmat(sol.o,m,1);
        % This is sqrt( -2*U'*V + a_i + b_j ) + o
        %ztmp = ztmp.*ok;
        
        %utmp = sqrt(abs(sol.w*sol.v+ repmat(sol.c_anchor,1,n)+repmat(sol.d2_anchor,m,1)))+repmat(sol.o,m,1);
        
        zcalc = zeros(size(z));
        zcalc(sol.rows,sol.cols) = ztmp;
        
        zok = zeros(size(z));
        zok(sol.rows,sol.cols) = ok;
    case {'rso'}
        z = sol.z;
        
        tmp = toa_calc_d_from_xy(sol.r, sol.s);
        ztmp = tmp + repmat(sol.o,m,1);
        
        zcalc = zeros(size(z));
        zcalc(sol.rows,sol.cols) = ztmp;
        zok = sol.inlmatrix;
        
end

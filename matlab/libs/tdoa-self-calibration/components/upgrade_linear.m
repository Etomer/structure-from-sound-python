function solut = upgrade_linear(sol)
% UPGRADE_LINEAR Upgrade a relaxed solution
%   sol = UPGRADE_LINEAR(sol) upgrades a relaxed solution, a solution in u,
%       v, a, b and o, to a solution in r, s and o. This is done linearly
%       using all data.

    D1J = sol.b - sol.b(1);
%     DI1 = sol.c_anchor;
    R = sol.u';
    S = sol.v;

    %S(:,j)'*H*S(:,j)-2*b*S(:,j) = D1J

    Hbase{1} = [1 0 0;0 0 0;0 0 0];
    Hbase{2} = [0 1 0;1 0 0;0 0 0];
    Hbase{3} = [0 0 1;0 0 0;1 0 0];
    Hbase{4} = [0 0 0;0 1 0;0 0 0];
    Hbase{5} = [0 0 0;0 0 1;0 1 0];
    Hbase{6} = [0 0 0;0 0 0;0 0 1];
    bbase{1} = [1 0 0]';
    bbase{2} = [0 1 0]';
    bbase{3} = [0 0 1]';
    A = zeros(size(S,2),9);
    B = zeros(size(S,2),1);
    for j = 2:size(S,2)
        for kk = 1:6
            A(j,kk)= S(:,j)'*Hbase{kk}*S(:,j);
        end
        for kk = 1:3
            A(j,kk+6)= -2*bbase{kk}'*S(:,j);
        end
        B(j,1)=D1J(j);
    end
    xxx = A\B;
    H = zeros(3,3);
    b = zeros(3,1);
    for kk = 1:6
        H = H+Hbase{kk}*xxx(kk);
    end
    for kk = 1:3
        b = b+bbase{kk}*xxx(6+kk);
    end

    % Sometimes the matrix H is not positive 
    % definite like it should.
    % This is a hack-fix
    minl = min(eig(H));
    if minl < eps
        H = H+(1.1*(eps-minl))*eye(size(H));
    end
    L = chol(inv(H));

    %
    solut.type = 'rso';
    solut.rows = sol.rows;
    solut.cols = sol.cols;
    solut.inlmatrix = sol.inlmatrix;
    solut.z = sol.z;
    %
    solut.o = sol.o;
    solut.r = L*(R+repmat(b,1,size(R,2)));
    solut.s = L'\S;
end


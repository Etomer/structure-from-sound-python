function problem = createRandomUpgradeProblem(m,n,d,sigma,weighting)
% CREATERANDOMUPGRADEPROBLEM Create a random low-rank upgrade problem
%   problem = CREATERANDOMUPGRADEPROBLEM(m,n) creates a random low-rank
%       problem (U,V,a,b,c) with ground truth data (Rgt,Sgt) consisting of
%       m receivers and n senders in 3D.
%   problem = CREATERANDOMUPGRADEPROBLEM(m,n,d,sigma,weighting) provides
%       additional arguments for specifying the dimension of the space, d,
%       the standard deviation of the noise in the distance measurements,
%       sigma, and the weighting used when creating the double compaction
%       matrix.
    
    if nargin < 3
        d = 3;
    end
    if nargin < 4
        sigma = 0;
    end
    if nargin < 5
        weighting = 'canonical';
    end

    % Random points in unit cube.
    Rgt = rand(d,m);
    Sgt = rand(d,n);

    % Random points on unit sphere.
%     Rgt = randn(3,m);
%     Rgt = Rgt./vecnorm(Rgt);
%     Sgt = randn(3,n);
%     Sgt = Sgt./vecnorm(Sgt);

    % Make sphere ellipsoids.
%     Rgt = randn(3)*Rgt; 
%     Rgt = Rgt+randn(3,1);
%     Sgt = randn(3)*Sgt;
%     Sgt = Sgt+randn(3,1);

    % Paraboloid.
%     t = randn(2,m);
%     Rgt(1:2,:) = t;
%     Rgt(3,:) = sum(t.^2);
%     Sgt(1:2,:) = t;
%     Sgt(3,:) = sum(t.^2);
    
    % Hyperboloid.
%     t = 2*pi*rand(1,m);
%     s = randn(1,m);
%     hyp = randn(4,1);
%     Rgt = [hyp(1)*sqrt(s.^2+hyp(4).^2).*cos(t); hyp(2)*sqrt(s.^2+hyp(4).^2).*sin(t); hyp(3)*s];
%     Sgt = randn(3,n);

    % Plane.
%     Rgt = rand(3,m);
%     Rgt(3,:) = rand();
%     Sgt = rand(3,n);
    
    % Random points on random 2D manifold.
%     M = randn(3,6);
%     M(:,4:5) = 0; % Will make linear system rank deficient.
%     s = randn(1,m);
%     t = randn(1,m);
%     Rgt = M*[s.^2; s.*t; t.^2; s; t; ones(1,m)];
%     s = randn(1,n);
%     t = randn(1,n);
%     Sgt = M*[s.^2; s.*t; t.^2; s; t; ones(1,n)];

    % Create distance matrix.
    Dgt = pdist2(Rgt',Sgt');
    D = Dgt+sigma*randn(size(Dgt));
    D2 = D.^2;

    % Define weighting.
    switch weighting
        case 'canonical'
            wr = zeros(m,1);
            wr(1) = 1;
            ws = zeros(n,1);
            ws(1) = 1;
        case 'uniform'
            wr = ones(m,1)/m;
            ws = ones(n,1)/n;
        case 'random'
            wr = randn(m,1);
            wr = wr/sum(wr);
            ws = randn(n,1);
            ws = ws/sum(ws);
        case 'distanceInv'
            wr = 1./sum(D,2);
            wr = wr/sum(wr);
            ws = 1./sum(D,1)';
            ws = ws/sum(ws);
        case 'distance'
            wr = sum(D,2);
            wr = wr/sum(wr);
            ws = sum(D,1)';
            ws = ws/sum(ws);
        case 'distance2Inv'
            wr = 1./sum(D2,2);
            wr = wr/sum(wr);
            ws = 1./sum(D2,1)';
            ws = ws/sum(ws);
        case 'distance2'
            wr = sum(D2,2);
            wr = wr/sum(wr);
            ws = sum(D2,1)';
            ws = ws/sum(ws);
        otherwise
            error('''%s'' is not a valid weighting.',weighting);
    end

    a = wr'*D2;
    b = D2*ws;
    c = wr'*D2*ws;

    Cr = eye(m)-wr;
    Cs = eye(n)-ws;

    M = -Cr'*D2*Cs/2;
    [U,S,V] = svd(M);
    U = S(1:d,1:d)*U(:,1:d)';
    V = V(:,1:d)';
    
    % Add random transformation.
%     H = randn(3);
%     U = H'\U;
%     V = H*V;
    
    problem.U = U;
    problem.V = V;
    problem.a = a;
    problem.b = b;
    problem.c = c;
    problem.Rgt = Rgt;
    problem.Sgt = Sgt;
    problem.Dgt = Dgt;
    problem.Dmeas = D;
    problem.M = M;
end


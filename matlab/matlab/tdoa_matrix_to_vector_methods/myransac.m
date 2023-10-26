function [bestl,bestinl]=myransac(usel,ransac_k,ransac_tol);
%%

maxinl = 0;
N = size(usel,2);
bestl = NaN*ones(3,1);
bestinl = [];
if N>5,
    uu = [usel(1:2,:);ones(1,N)];
    for i = 1:ransac_k,
        isel = randperm(N,2);
        l = cross(uu(:,isel(:,1)),uu(:,isel(:,2)));
        l = l/l(1);
        %figure(2); plot(abs(l'*uu),'*')
        err = abs(l'*uu);
        inl = find(err<ransac_tol);
        if length(inl)>maxinl,
            %figure(2); plot(abs(l'*uu),'*')
            %length(inl)
            maxinl = length(inl);
            bestl = l;
            bestinl = inl;
        end
    end
else
    bestl = NaN*ones(3,1);
    bestinl = [];
end
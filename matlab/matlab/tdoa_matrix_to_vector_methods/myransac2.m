function [uext]=myransac2(uext,ransac_k,ransac_tol,width,jjump,inl_min,lutn_max);
%%

n = max(uext(2,:));
%width = 100;
whalf = round(width/2);
allinl = zeros(1,size(uext,2));
for midj = whalf:jjump:(n-whalf),
    jsel = find( abs(uext(2,:)-midj) < width/2 );
    usel = uext(:,jsel);
    [bestl,bestinl]=myransac(usel,ransac_k,ransac_tol);
    %[midj length(bestinl)]
    if isfinite(bestl(1,1)),
        if length(bestinl)>inl_min,
            if abs(bestl(2))<lutn_max,
                allinl(jsel(bestinl))=ones(1,length(bestinl));
            end;
        end
    end
end

uext = uext(:,find(allinl));

function toas = tdoamatrix2tdoavector(u);

m = size(u,1);
n = size(u{1,2},2);
toas = zeros(m,n);
refch = 6;

for k = setdiff(1:m,refch),   
    utmp = u{refch,k};
    n = size(utmp,2);
    ktmp = repmat((1:4)',1,n);
    jtmp = repmat(1:n,4,1);
    uext0 = [utmp(:) jtmp(:) ktmp(:)]';
    
    nollor = find((abs(uext0(1,:))<3));
    uext0(:,nollor)=[];

    [uext1]=myransac2(uext0,100,10,100,10,10,3);
    [uext2]=myransac2(uext1,100,5,100,10,10,3);
    [uext3]=myransac2(uext2,100,3,20,3,6,3);    
    
    for j = 1:n,
        tmp = find(uext3(2,:)==j);
        if length(tmp)==1,
            toas(k,j)=uext3(1,tmp);
        elseif length(tmp)>1,
            [minv,mini]=min(uext3(3,tmp));
            tmp = tmp(mini);
            toas(k,j)=uext3(1,tmp);
        else
            toas(k,j)=NaN;
        end
    end
    
end;

function asol_out = more_spoints_v5_bundle_all_channels(asol,u,speedofsound,a_sr,threshold_dist,model);

if nargin<5,
    threshold_dist = 0.1;
end;

threshold_dist

%keyboard;

%% First find data in u that are sufficiently close to the reprojected
% measurements u_calc.
threshold_nr = 0;
s_in = asol.s;
r = asol.r;
nrinls_out = zeros(1,size(s_in,2));

jj = 1:size(s_in,2);
jj = jj(find(isfinite(s_in(1,:))));

I1 = [];
I2 = [];
J = [];
K = [];
U = [];
UC = [];
nrinls_out = zeros(1,size(s_in,2));
%keyboard;
%
for j = 1:size(s_in,2),
    if round(j/1000)==(j/1000)
        fprintf('%s',['-' num2str(j/1000)]);
    end

    
    if 1, % Get the data from nearby time instances
        jjj = j;
        udata = [];
        for i1 = 1:11,
            for i2 = (i1+1):12,
                udata = [udata [i1*ones(1,4);i2*ones(1,4);u{i1,i2}(:,jjj)';(jjj-j)*ones(1,4)]];
            end
        end
    end;
    behall = find(isfinite(udata(3,:)));
    udata = udata(:,behall);
    behall = find(udata(3,:)~=0);
    udata = udata(:,behall);
    udata(3,:)=udata(3,:)*speedofsound/a_sr;
    
    %jk = 3134;
    etts = s_in(:,j);
    
    if isfinite(etts(1,1)),
        reproj = sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(2,:))).^2 )) - ...
            sqrt(sum( (repmat(etts,1,size(udata,2))-r(:,udata(1,:))).^2 ));
        res = reproj-udata(3,:);
        nrinl = sum( abs(res)<threshold_dist );
        if nrinl>= threshold_nr,
            inl = find(abs(res)<threshold_dist);
            udata = udata(:,inl);
            reproj = reproj(:,inl);
            nrinls_out(1,j)=nrinl;
            I1 = [I1 udata(1,:)];
            I2 = [I2 udata(2,:)];
            J  = [J j*ones(1,size(udata,2))];
            K = [K 1*ones(1,size(udata,2))];
            U = [U udata(3,:)];
            UC = [UC reproj];
        end
    end
end;


asol_out.nrinls = nrinls_out;
asol_out.r = asol.r;
asol_out.s = asol.s;
asol_out.I1 = I1;
asol_out.I2 = I2;
asol_out.J = J;
asol_out.K = K;
asol_out.U = U;
asol_out.UC = UC;


%% Omvandla från r s till x y

%ures = asol4b.U-asol4b.UC;

bsol4b.x = asol_out.r;
bsol4b.y = asol_out.s;
data4b = [asol_out.I1;asol_out.I2;asol_out.J;asol_out.U]';
debug = 0;
xinorm = [2 6 12];

%% % cut away the first points, where the sound source is not located
% jcut = 740
% % Ändra i Bsol och data så att alla s-punkter från 1-740 tas bort och att
% % numreringen börjar med 1 från 741. 
% % Ändra samtidigt i data så att index för J sänks med 740.
% 
% bsol4b.y = bsol4b.y(:,(jcut+1):end);
% data4b(3,:)=data4b(:,3)-jcut;
% ok = find(data4b(:,3)>0);
% data4b = data4b(:,ok);

%% Design the smoothness connections cc and remove points in y that are NaN

y = bsol4b.y;
ok = isfinite(y(1,:));
NN = size(y,2);
j1= (1:(NN-2))';
j2 = j1+1;
j3 = j2+1;
cc = [j1 j2 j3];
okcc = find(all(ok(cc),2));
cc = cc(okcc,:);

converter1 = find(ok);
converter2 = zeros(1,NN);
converter2(converter1)=1:length(converter1);

% Remove y elements that are NaN. Remember to also change indices 
% in cc and
% in data4b(:,3);

cc2 = converter2(cc);
data4b(:,3)=converter2(data4b(:,3));
bsol4b.y = bsol4b.y(:,converter1);

%%
[bsol4bopt,res,jac]=bundletdoasmooth_allchannels(data4b,bsol4b,model,debug,xinorm,cc2);

%% Put the data back as r,s and also expand to the NaN elements of y

s_out = asol.s;
s_out(:,converter1)=bsol4bopt.y;

asol_out.r = bsol4bopt.x;
asol_out.s = s_out;
asol_out.res = res;
asol_out.jac = jac;


%% Demo script that reads one experiment with 12 recordings and runs the
% full system
addpath(genpath(".."))
%addpath("../")
%addpath("tdoa_matrix_to_vector_methods/")
%addpath("robust_multilateration_methods/")
speedofsound = 343.2;

%addpath("/Users/eriktegler/phd_studies/Forskning/tdoa_datasets_structuring/src/matlab/")

%% Set folder names


experiment_path = "../data/tdoa_20201016/data/music_0014/";
% ex.loadfiledir = fullfile(getenv('SFS_ROOT'), 'data', 'music');
% ex.savefiledir = fullfile(getenv('SFS_ROOT'), 'result');
% ex.XXXX = '0014';
% ex.loadname_raw   = fullfile(ex.loadfiledir,['music_' ex.XXXX '_raw_sound']);
% ex.loadname_gt   = fullfile(ex.loadfiledir,['music_' ex.XXXX '_gt']);
% ex.loadname_raw = strcat("../",ex.loadname_raw);

%% Read sound files

%load(ex.loadname_raw);
raw = read_experiment_audio(experiment_path);

%%


%% Correlation: GCC-PHAT. Find detections

settings.nbrOfSamples = size(raw.aaint,1);
settings.wf = @(x) 1./(abs(x)+(abs(x)<5e-3)); %weighting function

settings.firstSamplePoint = 1; %center sample point of first frame
settings.frameSize = 6048;     %width of frame in sample points
settings.dx = 1000;            %distance between frames in sample points
settings.frameOverlap = settings.frameSize-settings.dx; %overlap between frames
settings.sw = 800;             %clipping of search window
%Default: [@(x) 1./(abs(x)+(abs(x)<5e-3)),1,2048,1048,800]
settings.mm = 12;
settings.sr = raw.a_sr;
settings.channels = 1:12;



% folder = "../data/music_0008/";
% 
% for i=1:12
%     temp = audioread(strcat(folder,"Track ",num2str(i),".wav"));
%     if i == 1
%         raw.aaint = zeros(length(temp),12);
%     end
%     raw.aaint(:,i) = temp;
% end
%
%scores = gccscores_erik(double(raw.aaint)'/(2^24),settings);
scores = gccscores(double(raw.aaint)',settings);

%T = size(raw.aaint,1)/96000;
%tn = size(scores{1,2},2);
%tt = T*(0:(tn-1))/(tn-1);
%
temp = vector2frames(1:length(raw.aaint),settings.frameSize,...
        settings.frameOverlap,settings.firstSamplePoint);
tt = temp(round(settings.frameSize/2),:)/raw.a_sr;

%% Get the top 4 peaks in the GCC-phat scores

settings.nbrOfPeaks = 4;       %max number of peaks
settings.minPeakHeight = 0;%0.01; %min value of local maxima
%Default: [4,0.01]
%keyboard;
u = getdelays(scores,settings);
%result.matchings.u = u;

%% Go from tdoa-matrix representation to tdoa-vector. 

toas = tdoamatrix2tdoavector(u);
% convert from samples (tdoa) to meter (z)
z = speedofsound*toas/raw.a_sr;
detections.toas = toas;
detections.z = z;
detections.sr = raw.a_sr;
detections.v = speedofsound;

% save(savename,'detections');
%% Run solver using the tdoa vector measurements icassp 2020 and icassp 2021 papers

toas = detections.toas;
z = detections.z;
a_sr = detections.sr;
speedofsound = detections.v;

% Only use those sound events for which there are measurements to at least
% 5 microphones

ztmp = z;
okcols = find(sum(isfinite(z))>=5); % Select those columns (times) that have at least 5 measurements
ztmp = z(:,okcols);
ztmp = ztmp+0.1; % Hack. Somehow the solver breaks down if there are zeros in the ztmp matrix. Why? 
                 % Add 0.1

[r, s, o, sol] = tdoa(ztmp, 'display', 'iter', 'sigma', 0.01);

%% Visualize result after initial TDOA estimation

zcalc=tdoa_calc_u_from_xyo(r,s,o);

figure(20);
plot(sum( abs(zcalc-ztmp)<0.02 ),'*');

oks = find(sqrt(sum(s.^2))<4);
oks = find(sum( abs(zcalc-ztmp)<0.02 )>7);
oks = find(sum( abs(zcalc-ztmp)<0.02 )>-2);
figure(4);
hold off
plot3(r(1,:),r(2,:),r(3,:),'g*');
hold on
%plot3(s(1,oks),s(2,oks),s(3,oks),'b-');
plot3(s(1,oks),s(2,oks),s(3,oks),'b*');
axis([-10 10 -10 10 -10 10])

figure(21);
plot(okcols(oks),s(:,oks)','*');
ylim([-5 5])

%% Refine sound source positions - eusipco 2021 paper

% Choose one of several possible methods, see EUSIPCO 2021 paper
% asol_out = estimate_spoints_ChanHo(r,u,raw.speedofsound,raw.a_sr);
% asol_out = estimate_spoints_velasco(r,u,raw.speedofsound,raw.a_sr);
% asol_out = estimate_spoints_random_1_truncl1(r,u,raw.speedofsound,raw.a_sr);
% asol_out = estimate_spoints_random_1_l1(r,u,raw.speedofsound,raw.a_sr);
% asol_out = estimate_spoints_random_10_truncl1(r,u,raw.speedofsound,raw.a_sr);
% asol_out = estimate_spoints_nonlinearleastsquares(r,u,raw.speedofsound,raw.a_sr);
asol_out = estimate_spoints_ransac_truncl2(r,u,raw.speedofsound,raw.a_sr); % Proposed ransac+opt method

figure(21);
plot(asol_out.s','*');
ylim([-5 5])

%% Use the time continuity

asol3 = more_spoints_v3_use_interval_1(asol_out,u,raw.speedofsound,raw.a_sr);
asol3b = more_spoints_v3_use_interval_1(asol3,u,raw.speedofsound,raw.a_sr);
asol3c = more_spoints_v3_use_interval_1(asol3b,u,raw.speedofsound,raw.a_sr);
model.smeas = 0.1;
model.smotion = 0.01;
threshold_dist = 0.05;
asol4a = more_spoints_v5_bundle_all_channels(asol3c,u,raw.speedofsound,raw.a_sr,threshold_dist,model);
asol4b = more_spoints_v5_bundle_all_channels(asol4a,u,raw.speedofsound,raw.a_sr,threshold_dist,model);
asol4c = more_spoints_v5_bundle_all_channels(asol4b,u,raw.speedofsound,raw.a_sr,threshold_dist,model);
asol4c = more_spoints_v5_bundle_all_channels(asol4c,u,raw.speedofsound,raw.a_sr,threshold_dist,model);

%% Load ground truth

%load(ex.loadname_gt); % Ground Truth gt
%gt.rgt = rgt;
%gt.sgt_resamp = sgt_resamp;
%gt.sgt_mocap = s_gt;
%gt.tt_mocap = tt_mocap;

gt = read_experiment_gt_positions(experiment_path);
gt.tt = tt;
gt.tn = length(tt);
gt.sgt_resamp = interp1(gt.tt_mocap, gt.s_gt,tt)';

% It is perhaps only reasonable to check the positions when the music has
% started. OK is a vector that indicates if there is music at this time
filenr = 14;
musikstarter = [0 , 0 , 0 , 0 , 0 , 0 , 0 , 458699 , 366220 , 259188 , 313408 , 380186 , 0 , 746817 , 405888 ];
musikstart = musikstarter(filenr);
musikstart_toas = round(musikstart/1000);
OK = isfinite(gt.sgt_resamp(1,:));
OK(1:musikstart_toas)=zeros(1,musikstart_toas);

%% Visualize/evaluate

asol_out2 = calculate_IJKU(asol_out,u,raw.speedofsound,raw.a_sr);
brief_display_solution(asol_out2,gt,OK);
brief_display_solution(asol3,gt,OK);
brief_display_solution(asol3b,gt,OK);
brief_display_solution(asol3c,gt,OK);
brief_display_solution(asol4a,gt,OK);
brief_display_solution(asol4b,gt,OK);
brief_display_solution(asol4c,gt,OK);

figure(1); evalres = plot_solution(asol_out2,gt);
figure(2); evalres = plot_solution(asol3,gt);
figure(3); evalres = plot_solution(asol3b,gt);
figure(4); evalres = plot_solution(asol3c,gt);
figure(5); evalres = plot_solution(asol4a,gt);
figure(6); evalres = plot_solution(asol4c,gt);

%evaluate_solution(asol_out2,gt,OK);

%% Visualize final results

r = asol4c.r;
s = asol4c.s;

% zcalc=tdoa_calc_u_from_xyo(r,s,o);
% 
% figure(20);
% plot(sum( abs(zcalc-ztmp)<0.02 ),'*');

oks = find(sqrt(sum(s.^2))<4);
%oks = find(sum( abs(zcalc-ztmp)<0.02 )>7);
%oks = find(sum( abs(zcalc-ztmp)<0.02 )>-2);
figure(99);
hold off
plot3(r(1,:),r(2,:),r(3,:),'g*');
hold on
plot3(s(1,oks),s(2,oks),s(3,oks),'b-');


%indx = 1:100
%plot3(s(1,oks(indx)),s(2,oks(indx)),s(3,oks(indx)),'b*');
%axis([-10 10 -10 10 -10 10])
%%
figure(99)
hold off
plot3(r(1,:),r(2,:),r(3,:),'g*');
hold on
view(3)
xlim([min(min(r(1,:)),max(s(1,:))), max(max(r(1,:)),max(s(1,:)))])
ylim([min(min(r(2,:)),max(s(2,:))), max(max(r(2,:)),max(s(2,:)))])
zlim([min(min(r(3,:)),max(s(3,:))), max(max(r(3,:)),max(s(3,:)))])
grid on;
h = animatedline('MaximumNumPoints', 100);

% Force a 3D view
view(3);

for k = 1:size(s,2)
    if oks(k)
        addpoints(h,s(1,(k)),s(2,(k)),s(3,(k)));
    
        drawnow
    end
    pause(0.01)
end


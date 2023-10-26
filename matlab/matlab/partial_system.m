
%specified in function later
filepath = "./data/tdoa_20201016/data/music_0014/tdoa_vector_music_0014_to_matlab.csv";
target_location = "./data"

% will be called from the projects main directory so we add path to the
% matlab files
addpath(genpath("./matlab"))

%read data
ztmp = csvread(filepath);

% some data prepping
ztmp = ztmp(2:end,2:end) + 0.1;
ztmp = -ztmp';
%okcols = find(sum(isfinite(ztmp))>=5);

%solving tdoa problem
[r, s, o, sol] = tdoa(ztmp, 'display', 'iter', 'sigma', 0.02);
s = s - mean(r')';
r = r - mean(r')';

% store results
writematrix("sender_positions.csv", s);


%% Visualize result after initial TDOA estimation


zcalc=tdoa_calc_u_from_xyo(r,s,o);

figure(20);
plot(sum( abs(zcalc-ztmp)<0.05 ),'*');

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





%% Visualize final results
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


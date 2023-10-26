%% Initialization.
dataindex = 1; % Change 1-7 for different datasets.

%% Load data.
datanames = {'data1.mat', 'data2.mat', 'data3.mat', 'data4.mat', ...
    'data5.mat', 'data6.mat', 'data7.mat'};
load(fullfile('data', datanames{dataindex}));

room_limits = [min(rgt,[],2) max(rgt,[],2)]' + [-1; 1];
room_limits(1,3) = 0;

%% Plot measurements.
% fig = figure(10);
% plot(z', 'LineWidth', 2);
% xlim([1,size(z,2)]);
% xlabel('Sound event j');
% ylabel('TDOA [m]');
% ax = gca;
% ax.FontName = 'Times';
% ax.FontSize = 12;
% ax.Position(1) = 0.05;
% ax.Position(3) = 0.9;
% fig.Position(3) = 3*fig.Position(4);

%% Run system.
[r, s, o, sol] = tdoa(z, 'display', 'iter', 'sigma', 0.03, 'mul5', 3);

%% Register solution against ground truth.
% Scaling = false - use current speed measurement.
% Scaling = true - reestimate speed of sound.
[~, ~, transform] = procrustes(rgt', r', 'Scaling', true);
Q = transform.b * transform.T';
t = transform.c(1, :)';

r = Q * r + t;
s = Q * s + t;

fprintf('Speed of sound: %.2f m/s\n', sound_speed*transform.b);

%% Plot 3D.
figure(2);
plot3(rgt(1, :), rgt(2, :), rgt(3, :), 'k*');
hold on
plot3(sgt(1, :), sgt(2, :), sgt(3, :), 'k-');
plot3(r(1, :), r(2, :), r(3, :), 'go', 'Color', lines(1));
plot3(s(1, :), s(2, :), s(3, :), 'g-', 'Color', lines(1));
axis equal

cube(room_limits, 'k:');
hold off

xlim(room_limits(:, 1)');
ylim(room_limits(:, 2)');
zlim(room_limits(:, 3)');

legend('GT r', 'GT s', 'Est. r', 'Est. s', 'Location', 'NorthEast');
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 12);

%% Plot residuals.
plot_tdoa_res(sol);

%% Plot 2D overhead view.
figure(3);
plot(rgt(2, :), -rgt(1, :), 'k*');
hold on
plot(sgt(2, :), -sgt(1, :), 'k-');
plot(r(2, :), -r(1, :), 'go', 'Color', lines(1));
plot(s(2, :), -s(1, :), 'g-', 'Color', lines(1));
axis equal
hold off

xlim([-4 4]);
ylim([-3 2]);
% xlim(room_limits(:, 2)');
% ylim(room_limits(:, 1)');

legend('GT r', 'GT s', 'Est. r', 'Est. s', 'Location', 'NorthEast');
title('Top view');
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 12);

%% Print results.
fprintf('Microphone rms error: %f m\n', rms(vecnorm(rgt-r)));
fprintf('Missing data: %.0f %%\n', 100*nnz(~isfinite(z))/numel(z));
fprintf('Est. outlier: %.0f %%\n', 100*nnz(~sol.inlmatrix & isfinite(z))/numel(z));

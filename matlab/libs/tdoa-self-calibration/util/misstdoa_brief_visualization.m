function stats = misstdoa_brief_visualization(sol)
% MISSTDOA_BRIEF_VISUALIZATION Plot status report for current solution

[zcalc,zok] = misstdoa_calc_z(sol);
zerr = zcalc(:)-sol.z(:);
zinl = sol.inlmatrix(:);
zerr = zerr(find(zinl));

stats = [...
length(sol.rows) size(sol.z,1); ...
length(sol.cols) size(sol.z,2); ...
sum(sol.inlmatrix(:)) prod(size(sol.z)); ...
sum(zok(:)) prod(size(sol.z)); ...
std(zerr) 0];

detin = sol.inlmatrix;
detmiss = isnan(sol.z);
plot_im = ...
    (detmiss)*5 + ...       % missing data (black)
    (detin)*1 + ...         % true pos (green) - good
    (~detmiss & ~detin)*3;  % + ... % false pos (red) really bad
%    (truein & ~detin)*2 + ... % false neg (orange) quite bad
%    (trueout & ~detin)*4; % true neg (yellow) - good
plot_colmap = [0 1 0;0.2 0.7 0.2; 1 0 0;0.9 1 0; 0 0 0];
%
figure(1); clf; subplot(2,1,1);
image(plot_im);
colormap(plot_colmap);
title(['R:' num2str(length(sol.rows)) '/' num2str(size(sol.z,1)) ...
    ' C:' num2str(length(sol.cols)) '/' num2str(size(sol.z,2)) ...
    ' I:' num2str(sum(sol.inlmatrix(:))) '/' num2str(prod(size(sol.z))) ...
    ' P:' num2str(sum(zok(:))) '/' num2str(prod(size(sol.z)))  ...
    ' E: ' num2str(std(zerr))]);
subplot(2,1,2);
hist(zerr,100);
m = length(sol.rows);
n = length(sol.cols);
title(['Theta2 - Std: ' num2str(std(zerr)*length(zerr)/(length(zerr)-(3*m+3*n-6)))]);



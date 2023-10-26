function plot_tdoa_res(sol)

z = sol.z(sol.rows, sol.cols);
inlmatrix = sol.inlmatrix(sol.rows, sol.cols);

d = pdist2(sol.r', sol.s');
res = z - sol.o - d;

finiteres = res(isfinite(res));
inlres = res(inlmatrix);

figure(1001);
subplot(2, 1, 1);
histogram(finiteres);
title(sprintf('All residuals - mu=%.3f, std=%.3f', mean(finiteres), std(finiteres)));

subplot(2, 1, 2);
histogram(inlres);
title(sprintf('Inlier residuals - mu=%.3f, std=%.3f', mean(inlres), std(inlres)));

figure(1002);
inl = isfinite(res) & inlmatrix;
outl = isfinite(res) & ~inlmatrix;

[sortabsres, sortind] = sort(abs(res(:)));
inls = inl(sortind);
outls = outl(sortind);

plot(find(inls), sortabsres(inls), 'g.');
hold on
plot(find(outls), sortabsres(outls), 'r.');
hold off
title(sprintf('Sorted residuals (inl: %.1f %%, outl: %.1f %%, miss: %.1f %%)', ...
    100 * nnz(inl) / numel(inl), 100 * nnz(outl) / numel(outl), 100 * nnz(~isfinite(res)) / numel(res)));
legend('Inliers', 'Outliers', 'Location', 'NorthWest');

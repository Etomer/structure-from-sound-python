%% Generate synthetic data.
m = 15;
n = 30;
dim = 3;
sigma = 1e-2;
miss_ratio = 0.02;
out_ratio = 0.02;

[z, gt] = generate_synthetic_tdoa(m, n, dim, sigma, miss_ratio, out_ratio);

%% Run system.
[r, s, o, sol] = tdoa(z, 'display', 'iter', 'sigma', max(sigma, 1e-6));
% [r, s, o, sol] = tdoa_random(z, 'display', 'off', 'sigma', max(sigma, 1e-6), 'inits', 10);

%% Register result against ground truth.
[Q, t] = kabsch([r, s], [gt.r, gt.s]);

r = Q * r + t;
s = Q * s + t;

%% Plot results.
figure(2);
plot3(gt.r(1, :), gt.r(2, :), gt.r(3, :), 'r*');
hold on
plot3(gt.s(1, :), gt.s(2, :), gt.s(3, :), 'g*');
plot3(r(1, :), r(2, :), r(3, :), 'ro');
plot3(s(1, :), s(2, :), s(3, :), 'go');
axis equal
hold off

%% Print results.
er = vecnorm(r-gt.r);
es = vecnorm(s-gt.s);
eo = abs(o-gt.o);

fprintf('RMSE - r: %e, s: %e, o: %e\n', rms(er), rms(es), rms(eo));
fprintf('MAXE - r: %e, s: %e, o: %e\n', max(er), max(es), max(eo));

% Visualize the numerical stability of all solvers by plotting the
% residuals when providing noiseless data.

%% Initialization.
iters = 100;

%% Define solvers.
solvers = get_offset_solvers();

%% Run tests.
% Disable warnings when measuring execution time. This is restored below.
w = warning('off','all');

ress = nan(length(solvers), iters);
exectime = zeros(1, length(solvers));
for isolv = 1:length(solvers)
    solver = solvers(isolv);
    m = solver.m;
    n = solver.n;
    dim = solver.rank;

    fprintf('(%d/%d) Evalutating %s...\n', ...
        isolv, length(solvers), solver.name);

    totaltime = 0;
    for iter = 1:iters
        % Setup problem.
        r = randn(dim, m);
        s = randn(dim, n);
        o = randn(1, n);

        d = pdist2(r', s');
        z = d + o;

        % Run solver and measure execution time.
        tstart = tic();
        oest = solver.solv(z);
        totaltime = totaltime + toc(tstart);

        % Remove complex solutions.
        oest = real(oest(:, all(abs(imag(oest)) < 1e-6)));

        if isempty(oest)
            continue;
        end

        ress(isolv, iter) = min(vecnorm(oest - o', 2, 1));
    end

    exectime(isolv) = totaltime / iters;
end

% Restore warnings.
warning(w);

%% Plot numerical stability.
figure(1);
edges = -16:0.25:0;
centers = edges(1:end-1) + diff(edges) / 2;
for isolv = 1:length(solvers)
    solver = solvers(isolv);

    hc = histcounts(log10(ress(isolv, :)), edges) / iters;
    plot(centers, hc, 'LineWidth', 2);
    hold on

    fprintf('%s failed to find a real solution in %.1f %% of all tests.\n', ...
        solver.name, 100*sum(isnan(ress(isolv, :)))/iters);
end
hold off

legend(solvers.name, 'Interpreter', 'none');
xlabel('log_{10}(error)');

set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 12);

%% Plot execution time.
figure(2);
bp = bar(1000*exectime, 'FaceColor', 'flat');
bp.CData = lines(length(solvers));
ylabel('Execution time [ms]');
xticklabels({solvers.name});
xtickangle(45);
set(gca, 'TickLabelInterpreter', 'none');
title('Execution time');

fprintf('Execution times (microseconds):\n');
fprintf('%5.0f\n', 1e6*exectime);


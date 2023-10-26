% Visualize the numerical stability of all solvers by plotting the
% residuals when providing noiseless data.
% Second plot shows the number of times the solvers failed to find any
% real solution.

%% Setup.
m = 10;
n = 10;
iters = 100; % Increase for smoother histograms. 10000 in paper.
sigma = 0; % Noise in distance measurements.

solvers = getSolvers();

%% Solve random problems for all solvers.
ress = nan(length(solvers),iters);
failures = zeros(length(solvers),1);

% Define solver options (used for local optimization).
opts = struct();
opts.tol = 1e-9;
opts.maxIters = 20;
opts.refine = false;

for i=1:iters
    fprintf('Iteration: %d/%d\n',i,iters);
    
    % Contains the low rank matrices and ground truth measurement
    % matrix.
    prob = createRandomUpgradeProblem(m,n,3,sigma);
    
    for k=1:length(solvers)

        % Select solver.
        solver = solvers(k);

        % Create random minimal sample from the problem.
        sample = createRandomSample(solver.id,prob);
        
        % Solve minimal problem.
        [Lti,q] = solver.solve(sample,opts);

        if ~isempty(Lti)
            for j=1:length(Lti)
                Rhat = Lti{j}*sample.fullU;
                Shat = Lti{j}'\(sample.fullV+q(:,j));

                % Calculate error in distances.
                Dhat = pdist2(Rhat',Shat');
                ress(k,i) = min(ress(k,i),rms(prob.Dmeas(:)-Dhat(:)));
            end
        else
            failures(k) = failures(k)+1;
        end
    end
end

%% Plot residual distributions.
edges = -17:0.25:1; % Figure 1 in paper.
% edges = -16:0.25:2; % Figure 2 in paper.
centers = edges(1:end-1)+diff(edges)/2;

% Swap the 2nd and 3rd line color for consistency with Figure 3.
lineColors = lines(length(solvers));
tmp = lineColors(2,:);
lineColors(2,:) = lineColors(3,:);
lineColors(3,:) = tmp;

figure(1);
for k=1:length(solvers)
    hc = histcounts(log10(ress(k,:)),edges)/iters;
    if k <= 7
        plot(centers,hc,'-','LineWidth',2,'Color',lineColors(k,:));
    else
        plot(centers,hc,'--','LineWidth',2,'Color',lineColors(k,:));
    end
    hold on
end
hold off

title('Numerical stability for all solvers');
xlabel('log_{10}(error)');
% axis([edges(1) edges(end) 0 0.1]); % Figure 1 in paper.
axis([edges(1) edges(end) ylim]);
legend({solvers.name},'Location','EastOutside');

% set(gca,'FontName','Times');
% set(gca,'FontSize',14);

%% Plot failure rate.
figure(2);
bp = bar(failures/iters*100,'FaceColor','flat');
bp.CData = lineColors(1:length(solvers),:);
ylabel('%');
xticklabels({solvers.name});
title('Percentage of failures');


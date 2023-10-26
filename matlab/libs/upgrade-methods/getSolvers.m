function solvers = getSolvers()
% GETSOLVERS Get an array of all upgrade solvers
%   solvers = GETSOLVERS() gets an array of upgrade solvers that take a
%       minimal sample and returns the inverse transpose of L and q.
%
%       The name/id of the solver are, in order: 900, 810, 801, 720, 711,
%       630, 621, 540, 531, 450, 441.
%
%   See also createRandomSample.


    % These solvers do not require an action matrix solver.
    solvers(1).id = [9 0 0];
    solvers(1).solve = @upgrade900;
    solvers(2).id = [8 1 0];
    solvers(2).solve = @upgrade810;
    solvers(3).id = [8 0 1];
    solvers(3).solve = @upgrade801;
    
    % These are the same regardless of saturation.
    solvers(4).id = [7 2 0];
    solvers(4).solve = @(sample,opts) upgradeGeneral(sample,opts,@solver_toa_upgrade_720);
    solvers(5).id = [7 1 1];
    solvers(5).solve = @(sample,opts) upgradeGeneral(sample,opts,@solver_toa_upgrade_711);

    % These solvers have a finite number of solutions even when not
    % saturated. However, it seems like saturating results in a slightly
    % more numerically stable solver.
%     solvers(6).id = [6 3 0];
%     solvers(6).solve = @(sample,opts) upgradeGeneral(sample,opts,@solver_toa_upgrade_630);
%     solvers(7).id = [6 2 1];
%     solvers(7).solve = @(sample,opts) upgradeGeneral(sample,opts,@solver_toa_upgrade_621);
    solvers(6).id = [6 3 0];
    solvers(6).solve = @(sample,opts) upgradeGeneral(sample,opts,@solver_toa_upgrade_630_sat);
    solvers(7).id = [6 2 1];
    solvers(7).solve = @(sample,opts) upgradeGeneral(sample,opts,@solver_toa_upgrade_621_sat);
    
    % These have to be saturated to get a finite number of solutions.
    solvers(8).id = [5 4 0];
    solvers(8).solve = @(sample,opts) upgradeGeneral(sample,opts,@solver_toa_upgrade_540_sat);
    solvers(9).id = [5 3 1];
    solvers(9).solve = @(sample,opts) upgradeGeneral(sample,opts,@solver_toa_upgrade_531_sat);
    % We ignore the 450 solver and also solvers more nonlinear than 441 as
    % by swapping the receivers and senders one of the solvers here should
    % work and will probably be more stable.
%     solvers(10).id = [4 5 0];
%     solvers(10).solve = @(sample,opts) upgradeGeneral(sample,opts,@solver_toa_upgrade_450_sat);
    solvers(10).id = [4 4 1];
    solvers(10).solve = @(sample,opts) upgradeGeneral(sample,opts,@solver_toa_upgrade_441_sat);

    for k=1:length(solvers)
        solvers(k).name = sprintf('%d%d%d',solvers(k).id);
    end
end


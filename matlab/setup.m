%% Add paths
addpath(genpath('matlab'))
addpath(genpath(fullfile('libs','upgrade-methods')))
addpath(genpath(fullfile('libs','tdoa-self-calibration')))
setenv('SFS_ROOT', pwd);

%% Mex files
current_path = pwd;
setenv('EIGEN_DIR', fullfile(current_path ,'libs','eigen'));
cd(fullfile('libs', 'upgrade-methods'));
run('buildMexSolvers.m');
cd(current_path);
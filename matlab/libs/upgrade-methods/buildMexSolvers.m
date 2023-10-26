% Builds MEX files from the C++ solvers. It takes a few minutes.
eigenDir = getenv('EIGEN_DIR');
if isempty(eigenDir)
    error('Please set environment variable EIGEN_DIR or provide path.');
end

cppSolvers = dir('solvers/*.cpp');

for i = 1:length(cppSolvers)
    [~,name] = fileparts(cppSolvers(i).name);
    mexFileName = ['solvers/' name '.' mexext];
    if exist(mexFileName,'file')
        fprintf('(%d/%d) ''%s'' already exists. Skipping.\n',i,...
            length(cppSolvers),mexFileName);
    else
        fprintf('(%d/%d) Building %s...\n',i,length(cppSolvers),...
            cppSolvers(i).name);
        fname = fullfile(cppSolvers(i).folder,cppSolvers(i).name);
        mex(['-I"' eigenDir '"'],fname,'-outdir','solvers');
    end
end
fprintf('Done building solvers.\n');

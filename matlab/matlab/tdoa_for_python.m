function tdoa_for_python(folder)
   file = strcat(folder,"tdoa_vectors_to_matlab.csv");

    % will be called from the projects main directory so we add path to the
    % matlab files
    addpath(genpath("./matlab"));
    
    %read data
    ztmp = csvread(file);
    
    % some data prepping
    ztmp = ztmp(2:end,2:end) + 0.1;
    ztmp = ztmp';
    %okcols = find(sum(isfinite(ztmp))>=5);
    
    %solving tdoa problem
    [r, s, o, sol] = tdoa(ztmp, 'display', 'none', 'sigma', 0.02);
    s = s - mean(r')';
    r = r - mean(r')';
    
    % store results
    writematrix(s,strcat(folder,"sender_positions.csv"));
    writematrix(r,strcat(folder,"receiver_positions.csv"));
    writematrix(o,strcat(folder,"offsets.csv"));
end
function tdoa_for_python(filepath)
    addpath(genpath(".."))
    
    ztmp = csvread(strcat(filepath));
    ztmp = ztmp(2:end,2:end) + 0.1;
    ztmp = -ztmp';
    okcols = find(sum(isfinite(ztmp))>=5);
    
    

end
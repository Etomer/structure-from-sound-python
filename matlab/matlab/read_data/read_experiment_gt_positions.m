function gt = read_experiment_gt_positions(experiment_path)
    fpath = strcat(experiment_path, 'info.json'); 
    fid = fopen(fpath); 
    rawfile = fread(fid,inf); 
    str = char(rawfile'); 
    fclose(fid); 
    val = jsondecode(str);
    
    table = readtable(strcat(experiment_path,"gt_positions.csv"));
    %setup
    dims = ["x","y","z"];
    
    % read mic positions
    gt.rgt = zeros(3,val.number_of_mics);
    for mic = 1:val.number_of_mics
        for dimi = 1:3
            temp = table.(strcat("mic", num2str(mic), "_",dims(dimi)));
            gt.rgt(dimi,mic) = median(temp(isfinite(temp)));
        end
    end
    
    % read speaker position
    for dimi = 1:3
        temp = table.(strcat("speaker", "_",dims(dimi)));
        if dimi == 1
            gt.s_gt = zeros(length(temp),3);
        end
        gt.s_gt(:,dimi) = temp;
    end
    
    %read time
    gt.tt_mocap = table.time';
end
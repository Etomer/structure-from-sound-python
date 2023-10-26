function raw = read_experiment_audio(experiment_path)
    fpath = strcat(experiment_path, 'info.json'); 
    fid = fopen(fpath); 
    rawfile = fread(fid,inf); 
    str = char(rawfile'); 
    fclose(fid); 
    val = jsondecode(str);
    
    raw.speedofsound = 343.2;
    raw.a_sr = val.sampling_rate; 
    
    for i = 1:val.number_of_mics
        [y,fs] = audioread(strcat(experiment_path, "Track ", num2str(i), ".wav"));
        if fs ~= raw.a_sr
            warning("sampling rate of file doesn't match, info samplning rate")
        end
        if i == 1
            raw.aaint = zeros(length(y), val.number_of_mics);
        end
        raw.aaint(:,i) = y;
    end
end
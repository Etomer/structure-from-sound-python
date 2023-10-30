import scipy as sp
from glob import glob
import json
import numpy as np
import pandas as pd
import time
import os

# default constant values
CHUNK_LENGTH = 5000 # samples
SPEED_OF_SOUND = 343 # m/s
 
#def plot_gcc_phat_image(sounds, fs, tdoa_gt, times_gt, chunk_length=5000)

def interpolateGtAtTimes(recording_folder, interpolation_times):
    position_gt, _, times_gt = read_tdoa_sound_ground_truth(recording_folder)
    speakerPositions = position_gt["speaker"]
    interpolatedGt = [np.interp(interpolation_times, times_gt, speakerPositions[i]) for i in range(3)]
    return interpolatedGt



def evaluate_tdoa_estimator_on_recording(tdoa_estimator, recording_folder, print_status=True, chunk_length=CHUNK_LENGTH):
    if print_status:
        print(recording_folder)
        prev_status_update = time.time()
    fs, sounds = read_tdoa_sound(recording_folder)
    position_gt, tdoa_gt, times_gt = read_tdoa_sound_ground_truth(
        recording_folder)

    problems = divide_into_tdoa_chunks(sounds, fs, tdoa_gt, times_gt, chunk_length=chunk_length)

    nmics = problems[0][0].shape[0]
    nproblems = len(problems)

    tdoa_chunk_gt = np.zeros((nmics, nmics, nproblems))
    tdoa_chunk_estimation = np.zeros((nmics, nmics, nproblems))
    for i, problem in enumerate(problems):
        if print_status:
            nowtime = time.time()
            if prev_status_update + 1 < nowtime:
                print(str(i) + "/" + str(len(problems)))
                prev_status_update = nowtime
        tdoa_chunk_estimation[:, :, i] = tdoa_estimator(problem)
        tdoa_chunk_gt[:, :, i] = problem[1]
    return tdoa_chunk_estimation, tdoa_chunk_gt

def gcc_phat(problem, keep_dim = False, high_freq_cutoff=300):
    c =  sp.fft.fft(problem[0])
    nmics = problem[0].shape[0]

    if keep_dim:
        result = np.empty((nmics,nmics,problem[0].shape[1]))
    else:
        result = np.empty((nmics,nmics))
    for mic1 in range(problem[0].shape[0]):
        for mic2 in range(mic1+1, problem[0].shape[0]):
            temp = c[mic1,:]*np.conj(c[mic2,:])
            temp /= np.abs(temp) + 1e-10
            temp[high_freq_cutoff:-high_freq_cutoff] = 0
            res = sp.fft.ifft(temp)

            if keep_dim:
                result[mic1,mic2] = np.real(res)
            else: 
                bi = np.argmax(res)
                result[mic1,mic2] = bi if bi < problem[0].shape[1]/2 else bi - problem[0].shape[1]
            result[mic2,mic1] = -result[mic1,mic2]

    return result


def divide_into_tdoa_chunks(sounds, fs , tdoa_gt, times_gt, chunk_length=CHUNK_LENGTH, keep_missing_gt_chunks=False, speed_of_sound = SPEED_OF_SOUND):
    """
    Divides sound into chunks of a given length, and matches with ground truth for tdoa

    Input : 
        sounds, numpy_array of size=(nmics, sound_length)
        fs, int; sampling frequency
        tdoa_gt, numpy_array of size=(nmics, nmics, gt_len); containing time-difference groudn truth between pair of mics (in the unit samples)
        times_gt, numpy_array of size=(gt_len); contains times of ground_truth sampling, t=0 is when sound recording started

    Output : 
        list((
            sounds, numpy_array of size=(nmics, chunk_length)
            tdoa_gt, numpy_array of size=(nmics, nmics); 
                i.e. this assumes the chunk_length is short enough so that tdoa can be considered constant in the chunk
            time, int
                time since start of the original recording
            ))
    """

    problems = []
    starts = np.arange(0,sounds.shape[-1],chunk_length)[:-1]

    for start in starts:

        sound = sounds[:,start:(start+chunk_length)]

        gt_sample_in_period = (times_gt < (start + chunk_length)/fs) & (times_gt > start/fs)
        if (gt_sample_in_period.sum() == 0) and not keep_missing_gt_chunks:
            continue 

        gt_index = np.where(gt_sample_in_period)[0][0]

        problems.append((sound,tdoa_gt[:,:,gt_index]*fs/speed_of_sound, times_gt[gt_index]))
    return problems



def get_data_paths(data_folder, with_ground_truth=False):
    """
    returns paths to (experiments, recordings), where experiments are collection of recordings

    experiments : list(str) of experiment paths
    recordings : dict(str -> list(str)) with experiment paths as keys and list of recording paths as values
    """
    experiments = glob(data_folder + "/*/")
    
    #re_experiments = [exp.split("/")[-1] for exp in experiments]
    recordings = {exp:[] for exp in experiments}
    for i,exp in enumerate(experiments):
        temp = glob(exp + "/data/*")
        for j in temp:
            recordings[exp].append(j)

    if with_ground_truth:
        gt_recordings = {}
        for exp in experiments:
            temp = []
            for recording in recordings[exp]:
                info = json.load(open(recording + "/info.json",'r'))
                if info["has_ground_truth"]:
                    temp.append(recording)
            gt_recordings[exp] = temp
    
        recordings = gt_recordings


    return (experiments, recordings)


def read_tdoa_sound(case_folder):
    """
    Reads and packages a recording into a tuple containing sample rate and a numpy tensor with the shape (number_of_recordings, length_of_recordings)    
    """
    #sound_files = glob(case_folder + "/*.wav")
    
    info = json.load(open(case_folder + "/info.json",'r'))
    n = info["number_of_mics"]

    sounds = []
    for i in range(n):
        fs,s = sp.io.wavfile.read(case_folder + "/Track " + str(i+1)+".wav")
        sounds.append(s)

    return_obj = np.zeros((n,len(sounds[0])))
    for i in range(n):
        return_obj[i] = sounds[i]

    return (fs,return_obj)
    
    
def read_tdoa_sound_ground_truth(case_folder):
    """
    Reads and packages the ground truth to a recording
    
    RETURNS
    tuple : (positions, tdoa, time),
        positions : tuple {"mics":size=(number_of_mics, 3, timesteps), "speaker": size=(3, timesteps)}
        
        tdoa : size=(number_of_mics, numer_of_mics, timesteps),
            with tdoa(i,j) = toa_j - toa_i and tdoa(i,i) = nan, (in the units meters)
            
        time : (timesteps),
            time since recording started, for each of the timesteps
    """
    df_tdoa = pd.read_csv(os.path.join(case_folder, "gt_tdoa.csv"))
    df_pos  = pd.read_csv(os.path.join(case_folder, "gt_positions.csv"))
    info = json.load(open(os.path.join(case_folder, "info.json"),'r'))
    n = info["number_of_mics"]
    t = df_tdoa.shape[0]

    positions = {"mics": np.zeros((n,3,t)), "speaker": np.zeros((3,t))}
    
    tdoa = np.empty((n,n,t))
    time = np.array(df_tdoa['time'])

    dims = ['x','y','z']
    for i in range(n):
        for j,d in enumerate(dims):
            positions["mics"][i,j,:] = df_pos['mic' + str(i+1) + "_" + d]
    for j,d in enumerate(dims):
        positions["speaker"][j,:] = df_pos['speaker' + "_" + d]

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            else:
                tdoa[i,j,:] = df_tdoa['mic' + str(i+1) + "-mic" + str(j+1)]
                tdoa[j,i,:] = -df_tdoa['mic' + str(i+1) + "-mic" + str(j+1)]
    return (positions, tdoa, time)


    
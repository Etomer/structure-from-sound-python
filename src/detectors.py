from glob import glob
import numpy as np
import scipy as sp
import json
import os

from src.cnn_model_v2 import cnn_model_v2
import torch
# constants

SPEED_OF_SOUND = 343

# Helper functions ----------------------------


def _read_input_folder(input_folder):

    n = len(glob(os.path.join(input_folder, "Track*.wav")))
    # info = json.load(open(input_folder + "/info.json",'r'))
    # n = info["number_of_mics"]

    sounds = []
    for i in range(n):
        fs, s = sp.io.wavfile.read(os.path.join(input_folder, "Track " + str(i+1)+".wav"))
        sounds.append(s)

    return_obj = np.zeros((n, len(sounds[0])))
    for i in range(n):
        return_obj[i] = sounds[i]

    return (fs, return_obj)


def _divide_into_tdoa_chunks(sounds, fs, chunk_length):
    """
    Divides sound into chunks of a given length, and matches with ground truth for tdoa

    Input : 
        sounds, numpy_array of size=(nmics, sound_length)
        fs, int; sampling frequency

    Output : 
        list((
            sounds, numpy_array of size=(nmics, chunk_length)
            time, float
                time since start of the original recording in seconds
            ))
    """

    problems = []
    starts = np.arange(0, sounds.shape[-1], chunk_length)[:-1]

    for start in starts:

        sound = sounds[:, start:(start+chunk_length)]
        time = (start + chunk_length/2)/fs
        problems.append((sound, time))
    return problems


def _gcc_phat(problem, highest_fft_component_to_throw):
    temp = problem[0] - np.mean(problem[0],axis=1,keepdims=True)
    padded_problem = np.concatenate([temp, np.zeros(temp.shape)], axis=1)
    c = sp.fft.fft(padded_problem)
    nmics = problem[0].shape[0]

    result = np.empty((nmics, nmics))
    result[:] = np.nan
    for mic1 in range(padded_problem.shape[0]):
        for mic2 in range(mic1+1, padded_problem.shape[0]):
            temp = c[mic1, :]*np.conj(c[mic2, :])
            temp /= (np.abs(temp) + 1e-10)
            temp[highest_fft_component_to_throw:-
                 highest_fft_component_to_throw] = 0
            res = sp.fft.ifft(temp)

            bi = np.argmax(res)
            result[mic2, mic1] = bi if bi < padded_problem.shape[1] / \
                2 else bi - padded_problem.shape[1]
            result[mic1, mic2] = -result[mic2, mic1]

    return result


# Different detectors ----------------------------
def gcc_phat_detector(input_folder, output_folder=None, window_length=10000, speed_of_movement=1):
    """
    GCC-phat algorithm
    config : 
      window_length : #samples to use as window length
      speed_of_movement : max_speed the speaker is moving at (approximately) in m/s, (set this value to non-zero to start throwing out high frequency components)
    """
    if output_folder == None:
        output_folder = input_folder

    fs, sounds = _read_input_folder(input_folder)
    problems = _divide_into_tdoa_chunks(sounds, fs, chunk_length=window_length)

    highest_fft_component_to_throw = int(SPEED_OF_SOUND/speed_of_movement) if speed_of_movement != 0 else 0

    result = np.empty((len(problems), sounds.shape[0], sounds.shape[0]))
    result[:] = np.nan

    for i, problem in enumerate(problems):
        result[i,:,:] = _gcc_phat(problem, highest_fft_component_to_throw)*SPEED_OF_SOUND/fs
    times = [problem[1] for problem in problems]

    np.save(os.path.join(output_folder,"detections.npy"), result)
    np.save(os.path.join(output_folder,"detection_times.npy"), times)

def learned_detector(input_folder, output_folder=None, window_length=10000):
    if window_length != 10000:
        raise Exception("Wrong Window length")
    if output_folder == None:
        output_folder = input_folder

    fs, sounds = _read_input_folder(input_folder)
    #problems = _divide_into_tdoa_chunks(sounds[:,::6], fs//6, chunk_length=window_length)
    
    problems = []
    start_sep = window_length-3
    for start in np.arange(0,sounds.shape[1]-6*window_length,start_sep):
        problems.append((sounds[:,start:start+window_length*6:6],(start + window_length/2)/fs))

    model = torch.load("models/ResNetFFT/checkpoints/sonnet.pth", map_location="cpu", weights_only=True)

    result = np.empty((len(problems), sounds.shape[0], sounds.shape[0]))
    result[:] = np.nan

    for i, problem in enumerate(problems):
        
        #c = sp.fft.fft(problem[0])
        for mic_1 in range(result.shape[1]):
            for mic_2 in range(mic_1+1,result.shape[1]):
                #X = torch.tensor(np.stack([c[mic_1,:2500].real,c[mic_2,:2500].real,c[mic_1,:2500].imag,c[mic_2,:2500].imag],axis=0))
                X = torch.tensor(np.stack([problem[0][mic_1],problem[0][mic_2]],axis=0))
                X = X.unsqueeze(0).float()
                temp = 6*(model(X).argmax(dim=1) - 500)*SPEED_OF_SOUND/fs
                result[i,mic_1,mic_2] = -temp
                result[i,mic_2,mic_1] = temp
    times = [problem[1] for problem in problems]

    np.save(os.path.join(output_folder,"detections.npy"), result)
    np.save(os.path.join(output_folder,"detection_times.npy"), times)




def chirp_detector(input_folder, chirp_gt_path, output_folder = None):
    """
    First version of Chirp-detector
    """
    if output_folder == None:
        output_folder = input_folder

    #chirp_gt_path = "data/tdoa_20201016/meta/reference sound/Chirp 2 Hz.wav"
    fs_gt, sound_gt = sp.io.wavfile.read(chirp_gt_path)

    # cutout a single chirp from ground_truth
    zeros = np.where(np.abs(sound_gt[0:10000]) <= 1e-5*np.mean(abs(sound_gt)))[0]
    chirp_end_index = -1
    for i in range(len(zeros) - 5):
        if np.all(zeros[i+1:i+5] - zeros[i:i+4] - 1 == 0):
            chirp_end_index = zeros[i]
            break
    if chirp_end_index == -1:
        raise "couldn't find end of chirp"
    chirp = sound_gt[0:chirp_end_index]

    # load data 
    fs, sounds = _read_input_folder(input_folder)


    # resample (data or gt depending on which has lower sampling rate)
    sampling_fraction = fs/fs_gt
    if sampling_fraction > 1:
        sounds = np.stack([np.interp(np.arange(0,len(sound),sampling_fraction), np.arange(0,len(sound)), sound) for sound in sounds])
    else:
        sampling_fraction = 1/sampling_fraction
        chirp = np.interp(np.arange(0,len(chirp),sampling_fraction), np.arange(0,len(chirp)), chirp)

    # Find chirps in data and cut out them + some padding

    sound = sounds[0] # choose one of the mics to find the time of the chirps approximately
    # rolling average computation is much quicker than convolve
    cumsum_sound = np.cumsum(np.abs(sound))
    rolling_average = np.concatenate([np.zeros(len(chirp) // 2), cumsum_sound[len(chirp):] - cumsum_sound[:len(cumsum_sound) - len(chirp)], np.zeros(len(chirp) // 2)])

    def find_maxes(c, win_max_search = 10000, edge_trim = 10000):

        starts = np.arange(0,len(c), win_max_search)

        maxes = []
        for start in starts:
            cand = start + np.argmax(c[start:(start + win_max_search)])
            if np.argmax(c[(cand - win_max_search//2) : (cand + win_max_search//2)]) + (cand - win_max_search//2) == cand:
                if cand > edge_trim and len(c)-edge_trim > cand:
                    maxes.append(cand) 
        return np.array(maxes)

    maxes = find_maxes(rolling_average, win_max_search = len(chirp), edge_trim=len(chirp))

    problems = np.stack([sounds[:,maxes[max_i] - len(chirp) : maxes[max_i] + len(chirp)]/np.max(sounds[:,maxes[max_i] - len(chirp) : maxes[max_i] + len(chirp)]) for max_i in range(len(maxes))])
    chirp_pad = np.concatenate([np.zeros(len(chirp) // 2), chirp, np.zeros(len(chirp) // 2)])/np.max(chirp)

    def cross_correlation(signal, chirp):
        
        signal_pad = np.concatenate([np.zeros(len(signal)//2), signal, np.zeros(len(signal)//2)])
        chirp_pad = np.concatenate([np.zeros(len(signal)//2), chirp, np.zeros(len(signal)//2)])

        signal_tilde = sp.fft.fft((signal_pad))
        chirp_tilde = sp.fft.fft(np.flip((chirp_pad)))
        c = signal_tilde*chirp_tilde
        #c /= np.abs(c) + 1e-10 # GCC-PHAT
        score = sp.fft.ifft(c)
    
        return score

    shifts = []
    for mic_i in range(problems.shape[1]):
        temp_shift = []
        for i in range(problems.shape[0]):
            problem = problems[i,mic_i]
            score = cross_correlation(problem, chirp_pad)
            shift = np.argmax(score)
            if shift > len(score) / 2:
                shift = shift - len(score)
            temp_shift.append(shift)
        shifts.append(temp_shift)
    shifts = np.array(shifts)


    tdoa_vector = -SPEED_OF_SOUND*(shifts - shifts[0])/np.min([fs_gt,fs])
    detection_times = maxes/np.min([fs_gt,fs])
    detections = np.expand_dims(tdoa_vector.T,2) - np.expand_dims(tdoa_vector.T,1)
    for i in range(detections.shape[2]):
        detections[:,i,i] = np.nan
    # TODO: Detctions on diagonal should be nan
    np.save(os.path.join(output_folder,"detections.npy"),  detections)
    np.save(os.path.join(output_folder,"detection_times.npy"), detection_times)

    #np.save(os.path.join(output_folder,"tdoa_vectors.npy"), tdoa_vector)
    #np.save(os.path.join(output_folder,"tdoa_vector_times.npy"), detection_times)

    

    





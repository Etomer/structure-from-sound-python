from numpy import pi, sin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from glob import glob
import numpy as np
import scipy as sp
import json
import os,sys
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

if os.path.split(os.getcwd())[-1] != 'structure-from-sound-python':
    os.chdir("../")

sys.path.append('.')

from src.detectors import gcc_phat_detector
from src.tdoa_matrix_to_tdoa_vector import tdoa_matrix_to_tdoa_vector
import src.detectors as detectors


chirp_gt_path = "data/tdoa_20201016/meta/reference sound/Chirp 2 Hz.wav"
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
input_folder = "./data/tdoa_20201016/data/chirp_0002"
fs, sounds = detectors._read_input_folder(input_folder)

# resample (data or gt depending on which has lower sampling rate)
sampling_fraction = fs/fs_gt
if sampling_fraction > 1:
    sounds = np.stack([np.interp(np.arange(0,len(sound),sampling_fraction), np.arange(0,len(sound)), sound) for sound in sounds])
else:
    sampling_fraction = 1/sampling_fraction
    chirp = np.interp(np.arange(0,len(chirp),sampling_fraction), np.arange(0,len(chirp)), chirp)


sound = sounds[0] # choose one of the mics to find the time of the chirps approximately
# rolling average computation is much quicker than convolve
cumsum_sound = np.cumsum(np.abs(sound))
rolling_average = np.concatenate([np.zeros(len(chirp) // 2), cumsum_sound[len(chirp):] - cumsum_sound[:len(cumsum_sound) - len(chirp)], np.zeros(len(chirp) // 2)])

def find_maxes(c, win_max_search = 10000, edge_trim = 10000):

    starts = np.arange(edge_trim,len(c) - edge_trim, win_max_search)

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


i = 30
mic = 4
view_ind = np.arange(3100,9000)



fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, bottom=0.4)
axis_color = 'lightgoldenrodyellow'

shift_0 = 0
amp_0 = 0.2
shift_ref_0 = 0
amp_ref_0 = 0
amp_ref2_0 = 0
shift_ref2_0 = 0


shift_0_slider_ax  = fig.add_axes([0.15, 0.3, 0.65, 0.03], facecolor=axis_color)
shift_0_slider = Slider(shift_0_slider_ax, 'shift', -1200, 1200, valinit=shift_0)

amp_0_slider_ax  = fig.add_axes([0.15, 0.25, 0.65, 0.03], facecolor=axis_color)
amp_0_slider = Slider(amp_0_slider_ax, 'Amp', 0.0, 1.5, valinit=amp_0)

shift_ref_0_slider_ax  = fig.add_axes([0.15, 0.2, 0.65, 0.03], facecolor=axis_color)
shift_ref_0_slider = Slider(shift_ref_0_slider_ax, 'shift_ref', 0, 800.0, valinit=shift_ref_0)

amp_ref_0_slider_ax  = fig.add_axes([0.15, 0.15, 0.65, 0.03], facecolor=axis_color)
amp_ref_0_slider = Slider(amp_ref_0_slider_ax, 'Amp_ref', 0.0, 1.5, valinit=amp_ref_0)

shift_ref2_0_slider_ax  = fig.add_axes([0.15, 0.1, 0.65, 0.03], facecolor=axis_color)
shift_ref2_0_slider = Slider(shift_ref2_0_slider_ax, 'shift_ref', 0, 1200.0, valinit=shift_ref2_0)

amp_ref2_0_slider_ax  = fig.add_axes([0.15, 0.05, 0.65, 0.03], facecolor=axis_color)
amp_ref2_0_slider = Slider(amp_ref2_0_slider_ax, 'Amp_ref', 0.0, 1.5, valinit=amp_ref2_0)


sound = problems[i,mic,view_ind]
sound -= np.mean(sound)

def sliders_on_changed(val):
    ax.cla()

    ax.plot(sound)

    wav_dir = amp_0_slider.val*chirp_pad[view_ind - int(shift_0_slider.val)]    
    wav_ref = amp_ref_0_slider.val*chirp_pad[view_ind - int(shift_0_slider.val + shift_ref_0_slider.val)]
    wav_ref2 = amp_ref2_0_slider.val*chirp_pad[view_ind - int(shift_0_slider.val + shift_ref2_0_slider.val)]

    #ax.plot(wav_dir, 'b', alpha = 0.2)
    #ax.plot(wav_ref, 'r',alpha = 0.2)

    ax.plot(wav_dir + wav_ref + wav_ref2)
    fig.canvas.draw_idle()

shift_0_slider.on_changed(sliders_on_changed)
amp_0_slider.on_changed(sliders_on_changed)
shift_ref_0_slider.on_changed(sliders_on_changed)
amp_ref_0_slider.on_changed(sliders_on_changed)
shift_ref2_0_slider.on_changed(sliders_on_changed)
amp_ref2_0_slider.on_changed(sliders_on_changed)

plt.show()



# plt.figure(figsize=(16,6))
# plt.plot(gt_part)
# plt.plot(sound)




# def signal(amp, freq):
#     return amp * sin(2 * pi * freq * t)





# # Adjust the subplots region to leave some space for the sliders and buttons
# fig.subplots_adjust(left=0.25, bottom=0.25)

# t = np.arange(0.0, 1.0, 0.001)
# amp_0 = 5
# freq_0 = 3

# # Draw the initial plot
# # The 'line' variable is used for modifying the line later
# [line] = ax.plot(t, signal(amp_0, freq_0), linewidth=2, color='red')
# ax.set_xlim([0, 1])
# ax.set_ylim([-10, 10])

# # Add two sliders for tweaking the parameters

# # Define an axes area and draw a slider in it
# amp_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axis_color)
# amp_slider = Slider(amp_slider_ax, 'Amp', 0.1, 10.0, valinit=amp_0)

# # Draw another slider
# freq_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axis_color)
# freq_slider = Slider(freq_slider_ax, 'Freq', 0.1, 30.0, valinit=freq_0)

# # Define an action for modifying the line when any slider's value changes


# # Add a button for resetting the parameters
# reset_button_ax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
# reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
# def reset_button_on_clicked(mouse_event):
#     freq_slider.reset()
#     amp_slider.reset()
# reset_button.on_clicked(reset_button_on_clicked)

# # Add a set of radio buttons for changing color
# color_radios_ax = fig.add_axes([0.025, 0.5, 0.15, 0.15], facecolor=axis_color)
# color_radios = RadioButtons(color_radios_ax, ('red', 'blue', 'green'), active=0)
# def color_radios_on_clicked(label):
#     line.set_color(label)
#     fig.canvas.draw_idle()
# color_radios.on_clicked(color_radios_on_clicked)


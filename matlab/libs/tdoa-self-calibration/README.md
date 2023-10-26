# Fast and Robust Stratified Self-Calibration Using Time-Difference-of-Arrival Measurements
This repository contains suplementary materials associated with the paper [Fast and Robust Stratified Self-Calibration Using Time-Difference-of-Arrival Measurements](https://ieeexplore.ieee.org/abstract/document/9414309).

## Abstract
In this paper we study the problem of estimating receiver and sender positions using time-difference-of-arrival measurements. For this, we use a stratified, two-tiered approach. In the first step the problem is converted to a low-rank matrix estimation problem. We present new, efficient solvers for the minimal problems of this low-rank problem. These solvers are used in a hypothesis and test manner to efficiently remove outliers and find an initial estimate which is used for the subsequent step. Once a promising solution is obtained for a sufficiently large subset of the receivers and senders, the solution can be extended to the remaining receivers and senders. These steps are then combined with robust local optimization using the initial inlier set and the initial estimate as a starting point. The proposed system is verified on both real and synthetic data.

## Getting started
1. Clone this repository.
2. The upgrade solvers from [https://github.com/martinkjlarsson/upgrade-methods](https://github.com/martinkjlarsson/upgrade-methods) are needed. Clone that repository and run `buildMexSolvers.m` as instructed in the README. Place `upgrade-methods` and `tdoa-self-calibration` in the same folder, or make sure `upgrade-methods` is added to the path in MATLAB.
3. In MATLAB, run `startup.m` to setup the path.
4. Tests scripts for evaluating the minimal offset solvers and running the full TDOA system on real and sythetic data can be found in the `tests` folder.

## Data
Seven TDOA datasets can be found in the `data` folder. The setup consisted of 12 omni-directional microphones (the T-bone MM-1) spanning a volume of 4.0 × 4.6 × 1.5 meters. A speaker was moved through the setup while emitting sound, and ground truth positions for the microphones and speaker positions where found using a Qualisys motion capture system. A chirp sound was played with regular (dataset 1-5) or irregular (dataset 6-7) intervals. The arrival times were found using cross-correlation between the recordings and the original chirp. There was no missing data. The temperature in the room was measured to be 20.1 °C which indicates a speed of sound of 343 m/s. Each MAT-file contains the following variables:
* `name` - name of the dataset.
* `rgt` - ground truth receiver/microphone positions in meters.
* `sgt` - ground truth sender/speaker positions in meters. These do not correspond to the true positions of the chirp sound events, but represent a continuous track sampled at 150 Hz.
* `sound_speed` - estimated speed of sound (343 m/s) based on the temperature (20.1 °C) during data gathering.
* `z` - TDOA measurements in meters. Divide by `sound_speed` to get TDOA measurement in seconds.

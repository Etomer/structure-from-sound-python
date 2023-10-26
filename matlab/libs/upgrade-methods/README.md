# Upgrade Methods for Stratified Sensor Network Self-calibration
This repository contains suplementary materials associated with the paper [Upgrade Methods for Stratified Sensor Network Self-calibration](https://ieeexplore.ieee.org/abstract/document/9054025).

## Abstract
Estimating receiver and sender positions is often solved using a stratified, two-tiered approach. In the first step the problem is converted to a low-rank matrix estimation problem. The second step can be seen as an affine upgrade. This affine upgrade is the focus of this paper. In the paper new efficient algorithms for solving for the upgrade parameters using minimal data are presented. It is also shown how to combine such solvers as initial estimates, either directly or after a hypothesis and test step, in optimization of likelihood. The system is verified on both real and synthetic data.

## Getting started
1. Run `buildMexSolvers.m` to build the C++ files in `solvers/`. This takes a few minutes and requires [Eigen](http://eigen.tuxfamily.org). Set the environment variable `EIGEN_DIR` or modify `buildMexSolvers.m` to provided a path to Eigen.
2. Run `plotNumericalStability.m` to produce histograms showing the numerical stability of the solvers (see Figures 1 and 2 in the paper). Modify `createRandomUpgradeProblem.m` to constrain receiver and sender positions to a sphere or some other manifold.

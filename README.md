# Simulation Study: Comparing Adaptive Lasso and SCAD with Preliminary Shrinkage Estimators

## Description

This repository contains R code for a simulation study that compares the performance of Adaptive Lasso (ALASSO) and Smoothly Clipped Absolute Deviation (SCAD) methods, 
both with and without a Preliminary Shrinkage Estimator (PSE). The simulation evaluates the methods in terms of estimating both parametric (beta coefficients) and 
nonparametric components in a regression model.

## Contents

- `simulation.R`: The main R script containing all functions and code to run the simulation.
- `README.md`: This file.

## Requirements

The following R packages are required to run the simulation:

- `glmnet`
- `ncvreg`
- `ggplot2`
- `gridExtra`
- `reshape2`

You can install them using:

```R
install.packages(c("glmnet", "ncvreg", "ggplot2", "gridExtra", "reshape2"))

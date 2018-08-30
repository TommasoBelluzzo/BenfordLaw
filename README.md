# Benford Law

This script represents a full-featured framework for Benford's Law conformity assessment. It can be used in order to perform the following tests proposed by Nigrini et al. (2012):
* the `Primary Tests`: First Digits Analysis, Second Digits Analysis, First-Two Digits
* the `Advanced Tests`: Third Digits Analysis, Second Order Analysis, Summation Analysis
* the `Associated Tests`: Last-Two Digits Analysis, Number Duplication Analysis, Distortion Factor Model
* the `Mantissae Analysis`
* the `Zipf's Law Analysis`

For each significant digit analysis, the following conformity indicators are provided:

* **Goodness-of-Fit Measures (14):**
  * Anderson-Darling Discrete (Choulakian, 1994)
  * Chebyshev Distance (Leemis, 2000)
  * Cramer-von Mises Discrete (Choulakian, 1994)
  * Euclidean Distance (Cho & Gaines, 2007)
  * Freedman's U2 (Freedman, 1981)
  * Freeman-Tukey T2 (Freeman & Tukey, 1950)
  * Hotelling's Joint Digits (Hotelling, 1931)
  * Judge-Schechter Mean Deviation (Judge & Schechter, 2009)
  * Kolmogorov-Smirnov (Kolomonorgov, 1933)
  * Kuiper (Kuiper, 1960)
  * Likelihood Ratio (Neyman & Pearson, 1933)
  * Pearson's X2 (Pearson, 1900)
  * Watson's U2 (Choulakian, 1994)
* **Mean Absolute Deviation**
* **Sum of Square Differences**
* **Z-Scores**

In addition, the script is capable of producing Benford's Law conforming random numbers.

## Requirements

The minimum Matlab version required is `R2014a`. In addition, the `Statistics and Machine Learning Toolbox` must be installed in order to properly execute the script.

## Usage

1. Create a properly structured database (see the paragraph below).
1. Edit the `run.m` script following your needs.
1. Execute the `run.m` script.

## Dataset

This script doesn't require a dataset. Numeric data can be 

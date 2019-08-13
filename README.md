# Benford Law

This script represents a full-featured framework for assessing Benford's Law conformity. It can be used in order to perform all the tests proposed by [Nigrini et al. (2012)](https://www.nigrini.com/):

* the `Primary Tests`: First Digits Analysis, Second Digits Analysis, First-Two Digits
* the `Advanced Tests`: Third Digits Analysis, Second Order Analysis, Summation Analysis
* the `Associated Tests`: Last-Two Digits Analysis, Number Duplication Analysis, Distortion Factor Model
* the `Mantissae Analysis`
* the `Zipf's Law Analysis`

For each significant digit analysis, the following conformity indicators are provided:

* **Goodness-of-Fit Measures (14):**
  * Anderson-Darling Discrete ([Choulakian et al., 1994](https://www.researchgate.net/publication/245812965_Cramer-von_Mises_statistics_for_discrete_distributions))
  * Chebyshev Distance ([Leemis et al., 2000](https://www.researchgate.net/publication/2483540_Survival_Distributions_Satisfying_Benford's_Law))
  * Cramer-von Mises Discrete ([Choulakian et al., 1994](https://www.researchgate.net/publication/245812965_Cramer-von_Mises_statistics_for_discrete_distributions))
  * Euclidean Distance ([Cho & Gaines, 2007](https://www.researchgate.net/publication/243102672_Breaking_the_Benford_Law))
  * Freedman's U2 (Freedman, 1981)
  * Freeman-Tukey T2 (Freeman & Tukey, 1950)
  * Hotelling's Joint Digits (Hotelling, 1931)
  * Judge-Schechter Mean Deviation ([Judge & Schechter, 2009](https://www.researchgate.net/publication/228313715_Detecting_Problems_in_Survey_Data_Using_Benford's_Law))
  * Kolmogorov-Smirnov (Kolomonorgov, 1933)
  * Kuiper (Kuiper, 1960)
  * Likelihood Ratio (Neyman & Pearson, 1933)
  * Pearson's X2 (Pearson, 1900)
  * Watson's U2 Discrete ([Choulakian et al., 1994](https://www.researchgate.net/publication/245812965_Cramer-von_Mises_statistics_for_discrete_distributions))
* **Mean Absolute Deviation** ([Nigrini et al., 2012](https://www.nigrini.com/))
* **Sum of Square Differences** ([Kossovsky, 2014](https://www.researchgate.net/publication/316552141_Benford's_law_Theory_the_general_law_of_relative_quantities_and_forensic_fraud_detection_applications))
* **Z-Scores** ([Nigrini et al., 2012](https://www.nigrini.com/))

## Requirements

The minimum Matlab version required is `R2014a`. In addition, the `Statistics and Machine Learning Toolbox` must be installed in order to properly execute the script.

## Dataset & Usage 

The framework doesn't require any specific dataset structure. Numeric data can be extracted from any source or produced using any existing methodology, but a minimum amount of 1000 elements (with at least 50 unique observations) is required in order to perform a coherent analysis.

The `run.m` script provides an example of how this framework can be used, but all the functions located in the `Scripts` folder can be executed in standalone computation processes. It is recommended to validate and preprocess the dataset using the `benford_data` function. The `benford_analyse` functions can be used in order to perform a full automatic analysis of the dataset and plot the results. The `benford_random` function is an additional tool that produces random numbers whose digits follow the Benford's Law distribution.

## Screenshots

![Second Digits Analysis 1](https://i.imgur.com/fL3fbgO.png)

![Second Digits Analysis 2](https://i.imgur.com/rpRonnV.png)

![First-Two Digits Analysis 1](https://i.imgur.com/FDWDGBj.png)

![First-Two Digits Analysis 2](https://i.imgur.com/MEkU0pm.png)

![Mantissae Analysis](https://i.imgur.com/x0L5tqV.png)

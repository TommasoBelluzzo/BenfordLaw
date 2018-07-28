# Benford Law

This script represents a full-featured framework for Benford's Law conformity assessment. It can be used in order to perform the following tests proposed by Nigrini et al. (2012):
* the `Primary Tests`: First Digits Analysis, Second Digits Analysis, First-Two Digits
* the `Advanced Tests`: Third Digits Analysis, Second Order Analysis, Summation Analysis
* the `Associated Tests`: Last-Two Digits Analysis, Number Duplication Analysis, Distortion Factor Model
* the `Mantissae Analysis`
* the `Zipf's Law Analysis`

For each significant digit analysis, the following conformity indicators are provided:
* `Z-scores`
* `Mean Absolute Deviation`
* `Sum of Square Differences`
* 14 Goodness-of-Fit measures:

%         - AD for the Anderson-Darling test (Choulakian, 1994)
%         - CV for the Cramer-von Mises test (Choulakian, 1994)
%         - DC for the Chebyshev Distance test (Leemis, 2000)
%         - DE for the Euclidean Distance test (Cho & Gaines, 2007)
%         - FR for the Freedman's U2 test (Freedman, 1981)
%         - G2 for the Likelihood Ratio test (Neyman & Pearson, 1933)
%         - J2 for the Joenssen's J2 test (Joenssen, 2014)
%         - JD for the Hotelling's Joint Digits test (Hotelling, 1931)
%         - JS for the Judge-Schechter Mean Deviation test (Judge & Schechter, 2009)
%         - KS for the Kolmogorov-Smirnov test (Kolomonorgov, 1933)
%         - KU for the Kuiper test (Kuiper, 1960)
%         - T2 for the Freeman-Tukey T2 test (Freeman & Tukey, 1950)
%         - U2 for the Watson's U2 test (Choulakian, 1994)
%         - X2 for the Pearson's X2 test (Pearson, 1900)


## Requirements

The minimum Matlab version required is `R2014a`. In addition, the `Statistics and Machine Learning Toolbox` must be installed in order to properly execute the script.

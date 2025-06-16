# stroke_EEG
Collaboration with Jin-Moo Lee and Asher Albertson, secondary data analysis of age related changes in the EEG spectral slope in people with stroke and healthy controls.

This repository contains the statistical analysis scripts and summary output that support the following manuscript:
Alberston, A.J., Landsness, E.C., Eisfelder, M., Young, B.M., Judge, B., Brier, M.R., Euler, M.J., Cramer, S.C., Lee, J-M., & Lohse, K.R. (2024). Stroke is associated with regional and age-specific changes in periodic and aperiodic cortical activity. bioRxiv. https://doi.org/10.1101/2024.11.07.622359

Unfortunately, due to the data use agreements of this secondary data analysis, we are not at liberty to share the raw data that would be required to make the analyses fully reproducible. Hopefully, however, the code and summary data contain here provide some additional transparency to make it easier for interested groups to see how the data were handled and analyzed. Please direct any questions to Keith Lohse, PhD, PStat at lohse@wustl.edu

In brief, the repository contains 5 scripts:
1. A data management script for merging the "participant level" data as they come out of our spec_param pipeline implemented in MATLAB: https://github.com/margareteisfelder/Automated-EEG-Cleaning-Pipeline/blob/main/README.md

2. A script for harmonizing the control datasets.

3. A script for imputing exponents and offsets in the control data using mulitple imputation.
   - Note that periodic components were not included.
   

4. A script for imputing exponents and offsets in the data from people with stroke using multiple imputaiton.
   - Note that periodic components were not included.
   

5. A statistical analysis script that contains all of the analyses conducted and reported in our paper.
   - Comparisons between EEG data for people with stroke and health controls.
     - Aperiodic Exponent
     - Number of Peaks of Peaks Found
     - Central Frequency of Peaks Found
     - Power at the Central Frequency
  - Subgroup analyses in people with stroke specifically.
     - Aperiodic Exponent (ipsilesional versus contralesional, midline electrodes excluded)
     - Number of Peaks of Peaks Found (ipsilesional versus contralesional, midline electrodes excluded)
     - Central Frequency of Peaks Found (ipsilesional versus contralesional, midline electrodes excluded)
     - Power at the Central Frequency (ipsilesional versus contralesional, midline electrodes excluded)
     - Association of the Exponent to the Box and Block Test



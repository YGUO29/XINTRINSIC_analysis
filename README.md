# XINTRINSIC_analysis

Analysis code used for wide-field imaging (intrinsic & calcium), for XINTRINSIC system developed in Wang Lab @JHU
(in progress...)

### Pre-process of trial-based data: 
- script_TrialBased.m: run first 3 sections in this script to combine repetitions, perform normalization (calculate deltaF/F0) and visualize data.
- ViewData.m: for visualize wide-field imaging movies, and get data matrix.

### Registration
(need NoRMCorre https://github.com/flatironinstitute/NoRMCorre)
- script_Registration.m: get motion metrics, perform rigid registration selected moving frames  

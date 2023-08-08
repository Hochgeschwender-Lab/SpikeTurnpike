# SpikeTiming
Computational routines for calculating decoding quality, precision, and reliability of spike timing.
## Overview
The spike timing analysis implements the time-domain white-noise (WN) analyses for spike trains, yielding the results described in Levi, Spivak, Sloin, Someck, Stark, 2022, Cell Reports.

## Simulated spike trains
To demonstrate the entire pipeline for simulated spike trains, use st_sim.m.

## Neuronal spike trains
To demonstrate the entire pipeline for real neuronal spike trains, use st_neuronal.m.

## Routines

### Wrappers
- st_sim.m
  - A wrapper, calls stAnalysis routine for simulated spikes trains
- st_neuronal.m
  - A wrapper, calls stAnalysis routine for neuronal spikes trains
- stAnalysis.m
  - Does time-domain WN analysis for spike trains

### Analysis
- sta.m
  - Compute sta from spike train and analog vector
- stacalc.m
  - Compute multi-sta from tagged spike trains and an analog vector
- stahat.m
  - Reconstruct analog signal
- stcodertypes.m
  - Classify reliability profiles into rate or temporal coders
- stmix.m
  - Shuffle/jitter spikes in a sparse array
- streconstruct.m
  - Reconstruct analog waveform from spike trains
- streconstructXval.m
  - Cross-validation wrapper for streconstruct.m
- streliability.m
  - Estimate reliability of single-cell spike trains
- stwiener.m
  - Construct a modified Wiener filter for spike train and analog signal

### Utilities
- bayescount.m
  - Compute the number of bins for bias calculation, using the Bayes counting procedure of Panzeri&Treves 96
- calc_bias_PT.m
  - Calculate an analytical estimate of MI bias from a joint probability table
- calc_fwhh.m
  - Determine full-width at half-height
- calc_pearson.m
  - Compute Pearson's correlation coeffiecient
- calc_spearman.m
  - Compute Spearman's correlation coeffiecient
- isisort.m
  - Sort spikes into ISI groups
- mixmat.m
  - Mix matrix elements
- mutinf.m
  - Mutual information from empirical distributions/counts
- my_xcorr.m
  - Compute the column-wise cross-correlation between real matrices
- myjet.m
  - Modified jet with extreme values in pure R,B.
- nangeomean.m
  - Geometric mean
- ParseArgPairs.m
  - Flexible argument assigning
- plot_raster.m
  - Raster display for spike trains
- resampleSig.m
  - Wrapper for resample.m

### C Utilities
These utility functions are operating system specific and require compiling in MATLAB (see section 'To run the code' below). 
- calc_cor.c
  - Calculate delays by enumerating on first spike train only
- zcs.c
  - Compute the Z-score of the columns of a given matrix

### Data
- WN10k.mat
  - WN pattern
  - WN Fs
- sim_spikes.mat 	
  - Spike trains of four simulated units: High precision rate coder, high precision temporal coder, low precision rate coder, and low precision temporal coder.
- neuronal_spikes.mat 	
  - Spike trains of 3 recoded units: DA PYR, DA INT, and IDA INT.

## To run the code
- Download all routines and data
- Compile the C utility functions using the 'mex' command in MATLAB:
  - To compile calc_cor.c, in MATLAB, write: mex('calc_cor.c')
  - To compile zcs.c, in MATLAB, write: mex('zcs.c')
- For simulated spike trains, in MATLAB, write: st_sim 
- For real neuronal spike trains, in MATLAB, write: st_neuronal

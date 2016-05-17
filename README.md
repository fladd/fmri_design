# fMRI design

## About

This is a work in progress. The goal is to create a Python pipeline to create and evaluate fMRI desgins (i.e. stimulus timings).

Currently, two following modules are available:

* **`make_random_timing.py` - Create a random fMRI design**

  This script is part of [AFNI](https://afni.nimh.nih.gov/afni/) and as such licensed under GPL (see also https://afni.nimh.nih.gov/pub/dist/src/README.copyright).
  
  Have a look at the documentation (https://afni.nimh.nih.gov/pub/dist/doc/program_help/make_random_timing.py.html).

* **`timings2bv.py` - Create BrainVoyager .prt and .sdm files from a set of timing files**

  One timing file represents one experimental condition and includes a single
  line with times in msec (float), separated by spaces (`e.g. 2.0 4.0 6.0`).
  
  Each timing file has the following naming convention: `*_<condition_name>.*` (e.g. `myexp_condition1.txt`).
  
  Timing files can, for instance result from the `make_random_timing.py` script.

## Example

We will create .prt and .sdm files for one random single-run TR-locked (TR=2.0 seconds) fMRI design of 420 seconds run length, with 2 conditions ('congruent' and 'incongruent'), 100 repetitions per conditions (each lasting 1 sec).

1. Run `make_random_timing.py`:

    ```
    python make_random_timing.py -num_stim 2                          \
                                 -stim_dur 1                          \
                                 -num_runs 1                          \
                                 -run_time 420                        \
                                 -num_reps 100                        \
                                 -prefix stimes.001                   \
                                 -pre_stim_rest 6                     \
                                 -post_stim_rest 12                   \
                                 -min_rest 0                          \
                                 -stim_labels "congruent incongruent" \
                                 -seed 1234567                        \
                                 -tr 2.0                              \
                                 -tr_locked                           \
                                 -show_timing_stats                   \
    ```
    
2. Run `timings2bv.py`:

    ```
    python timings2bv.py -d 1 stimes.001_* 420 2.0 
    ```
    
This results in the files `Stimulation.prt` and `Stimulation.sdm`, which can be loaded into BrainVoyager for further inspection.

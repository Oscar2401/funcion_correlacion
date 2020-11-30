# Correlation Function.

It are codes in C++ and python for efficiently computing 2-point and 3-point correlation functions.

## Introduction.

In the last decades, the Correlation Functions (CF) have been widely used in Cosmology. These CF provide a rich and abundant information about how matter is distributed in the Universe. In this work, we provide a set of codes written in C++ to calculate 2-point and 3-point correlation functions of galaxy distributions. Also, the codes were parallelized to be run on GPUs. 

## Dependencies.

To compile and run these codes is necessary: the gcc C compiler and the GSL libraries. Although it is parallelized for shared with OpenMP so it is necessary the compiler extension.

The codes was written in C++ but Python is used to plot and show the outputs of each code. For thise reason, matplotlib and numpy libraries are indispensable.



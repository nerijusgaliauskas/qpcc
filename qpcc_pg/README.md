October 19, 2015

1. ABOUT

These folders contain files that implement and help to execute an algorithm, devoted to find a GLOBAL minimizer of an optimization problem arising in multidimensional scaling with city-block distances. The algorithm is based on a parallel branch and bound method and is implemented using MPI.

2. INPUT DATA

A set of dissimilarity matrices are presented in the folder /input/data. The file /input/input defines names of the dissimilarity matrices that are used by the algorithm.

3. COMPILATION AND EXECUTION

In order to compile the program, LAPACK library (see http://www.netlib.org/lapack/) and an implementation of MPI (for example, MPICH) are required. Path to the LAPACK library has to be set in the file /program/run.sh. In order to start the program, script /program/run.sh has to be executed in the following form:

$ /program/run.sh N

where N is the desired number of processes. The program was tested on a computer with Intel® Core™ i7 4710HQ Processor by using from 1 to 8 processes. Thus, in order to run the program on a cluster of computers, a particular script has to be written.

4. OUTPUT

After the execution, a certain folder is constructed. It contains an output file for every dissimilarity matrix separately.
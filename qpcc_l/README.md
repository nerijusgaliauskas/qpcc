December 27, 2014

1. ABOUT

These folders contain files that implement and help to execute an algorithm, devoted to find a LOCAL minimizer of an optimization problem, arising in multidimensional scaling with city-block distances. The algorithm is described in the following article:

Roger Fletcher, Nerijus Galiauskas*, Julius Žilinskas*.
Quadratic Programming with Complementarity Constraints for Multidimensional Scaling with City-Block Distances
Baltic Journal of Modern Computing.
Vol. 2 (2014), No. 4, pp. 248-259.

2. INPUT DATA

A set of dissimilarity matrices are presented in folder /input/data. File /input/input defines names of the dissimilarity matrices that are used by the algorithm.

3. COMPILATION AND EXECUTION

In order to compile the program, LAPACK library is required (see http://www.netlib.org/lapack/). Path to the LAPACK library has to be set in the file /program/run.sh. In order to start the program, script /program/run.sh has to be executed.

4. OUTPUT

After the execution, a certain folder is constructed. It contains an output file for every dissimilarity matrix.

*The research was funded by a Grant (No. MIP-063/2012) from the Research Council of Lithuania.
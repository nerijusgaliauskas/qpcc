#!/bin/bash

# the main variables
blaslapack=/path/to/lapack

cd ..
data=$(pwd)/input/data
input=$(pwd)/input/input
source=$(pwd)/source
program=$(pwd)/program
program_name=qpcc_pglobal
output=$(pwd)/output
if [ ! -d $output ]; then
	mkdir $output
fi

# compiles the program
cd $source
mpifort qpcc_mod.f90 qpcc_pglobal_dr.f90 -L$blaslapack -llapack -lrefblas
mv a.out $program/$program_name
rm *.mod

# executes the program
while read line; do
	if [ ${line:0:1} != "#" ]; then
		echo $program_name $data/$line
		mpiexec -n $1 $program/$program_name < $data/$line > $output/$line.out
	fi
done < $input

# removes unnecessary files
rm $program/$program_name

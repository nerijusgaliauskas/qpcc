#!/bin/bash

# the main variables
blaslapack=/path/to/LAPACK

cd ..
data=$(pwd)/input/data
input=$(pwd)/input/input
source=$(pwd)/source
program=$(pwd)/program
program_name=qpcc_global
output=$(pwd)/output
if [ ! -d $output ]; then
	mkdir $output
fi

# compiles the program
cd $source
gfortran qpcc_mod.f90 qpcc_global_dr.f90 -L$blaslapack -llapack -lrefblas
mv a.out $program/$program_name
rm *.mod

# executes the program
while read line; do
	if [ ${line:0:1} != "#" ]; then
		echo $program_name $data/$line
		$program/$program_name < $data/$line > $output/$line.out
	fi
done < $input

# removes unnecessary files
rm $program/$program_name

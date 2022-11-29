#!/bin/bash

test_name=$1
root_dir=${2:-.}
data_dir=${root_dir}/${test_name}
work_dir=./${test_name}

echo "--------------------------------"
echo "runtest ${test_name}"
echo "--------------------------------"
echo "data_dir: ${data_dir}"
echo "work_dir: ${work_dir}"
echo

# copy data to work directory
if [ ! $work_dir -ef $data_dir ]; then
    echo "copy data from $work_dir to $data_dir"
    mkdir -p $work_dir
    cp -r ${data_dir}/*.def ${data_dir}/*.dat ${data_dir}/*.in ${data_dir}/*.toml ${data_dir}/output_ref ${work_dir}
fi

cd ${work_dir}

# run program
python3 ${root_dir}/src/qlms.py input.toml

# show results
if [ -f output/energy.dat ]; then
    echo "----- output -----"
    cat output/energy.dat
    echo "----- reference -----"
    cat output_ref/energy.dat
    echo
fi

# check results
python3 ${root_dir}/tests/compfile.py output_ref/energy.dat output/energy.dat

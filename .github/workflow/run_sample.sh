#!/bin/sh

if [ $# -lt 2 ]; then
  echo "Usage: sh $0 sample_directory path/to/uhfdry_out"
  exit 1
fi
sample_dir=$(cd $1 && pwd)
uhfdry_out=$(readlink -f $2)

nfails=0

cd ${sample_dir}/hubbard_chain/UHF
hwave ./input.toml
if [ $? -ne 0 ]; then
  echo hubbard_chain/UHF fails.
  nfails=`echo "$nfails + 1" | bc`
fi

cd ${sample_dir}/hubbard_chain/UHFk
$uhfdry_out stan.in && hwave ./input.toml && python3 ./output_band.py
if [ $? -ne 0 ]; then
  echo hubbard_chain/UHFk fails.
  nfails=`echo "$nfails + 1" | bc`
fi

cd ${sample_dir}/Hubbard_square/UHF
$uhfdry_out stan.in && hwave ./input.toml
if [ $? -ne 0 ]; then
  echo Hubbard_square/UHF fails.
  nfails=`echo "$nfails + 1" | bc`
fi

cd ${sample_dir}/Hubbard_square/UHFk
$uhfdry_out stan.in && hwave ./input.toml
if [ $? -ne 0 ]; then
  echo Hubbard_square/UHFk fails.
  nfails=`echo "$nfails + 1" | bc`
fi

cd ${sample_dir}/CDW_SDW
python3 ./run.py $uhfdry_out
if [ $? -ne 0 ]; then
  echo CDW_SDW fails.
  nfails=`echo "$nfails + 1" | bc`
fi

if [ $nfails -gt 0 ]; then
  exit 1
else
  exit 0
fi

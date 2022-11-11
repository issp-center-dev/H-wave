#!/bin/sh

if [ $# -lt 2 ]; then
  echo "Usage: sh $0 sample_directory path/to/uhfdry_out"
  exit 1
fi
sample_dir=$(cd $1 && pwd)
uhfdry_out=$(readlink -f $2)

nfails=0
failed_dirs=""

sample_name="hubbard_chain/UHF"
echo start $sample_name
cd ${sample_dir}/${sample_name}
hwave ./input.toml
if [ $? -ne 0 ]; then
  echo ${sample_name} failed.
  failed_dirs="${failed_dirs} ${sample_name}"
  nfails=`echo "$nfails + 1" | bc`
else
  echo ${sample_name} finished.
fi

sample_name="hubbard_chain/UHFk"
echo start $sample_name
cd ${sample_dir}/${sample_name}
$uhfdry_out stan.in && hwave ./input.toml && python3 ./output_band.py
if [ $? -ne 0 ]; then
  echo ${sample_name} failed.
  failed_dirs="${failed_dirs} ${sample_name}"
  nfails=`echo "$nfails + 1" | bc`
else
  echo ${sample_name} finished.
fi

sample_name="Hubbard_square/UHF"
echo start $sample_name
cd ${sample_dir}/${sample_name}
$uhfdry_out stan.in && hwave ./input.toml
if [ $? -ne 0 ]; then
  echo ${sample_name} failed.
  failed_dirs="${failed_dirs} ${sample_name}"
  nfails=`echo "$nfails + 1" | bc`
else
  echo ${sample_name} finished.
fi

sample_name="Hubbard_square/UHF"
echo start $sample_name
cd ${sample_dir}/${sample_name}
$uhfdry_out stan.in && hwave ./input.toml
if [ $? -ne 0 ]; then
  echo ${sample_name} failed.
  failed_dirs="${failed_dirs} ${sample_name}"
  nfails=`echo "$nfails + 1" | bc`
else
  echo ${sample_name} finished.
fi

sample_name="CDW_SDW"
echo start $sample_name
cd ${sample_dir}/${sample_name}
python3 ./run.py $uhfdry_out
if [ $? -ne 0 ]; then
  echo ${sample_name} failed.
  failed_dirs="${failed_dirs} ${sample_name}"
  nfails=`echo "$nfails + 1" | bc`
else
  echo ${sample_name} finished.
fi

if [ $nfails -gt 0 ]; then
  for failed_dir in $failed_dirs; do
    echo $failed_dir
  done
  exit 1
else
  exit 0
fi

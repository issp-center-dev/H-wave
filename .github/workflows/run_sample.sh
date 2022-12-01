#!/bin/sh

if [ $# -lt 2 ]; then
  echo "Usage: sh $0 sample_directory path/to/uhfdry_out"
  exit 1
fi
sample_dir=$(cd $1 && pwd)
uhfdry_out=$(readlink -f $2)

nfails=0
failed_dirs=""

sample_name="Hubbard_chain/UHF"
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

sample_name="Hubbard_chain/UHFk"
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

sample_name="Hubbard_cubic/UHF"
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

sample_name="Hubbard_cubic/UHFk"
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


# sample_name="Hubbard_cubic/finiteT/UHFk"
# echo start $sample_name
# cd ${sample_dir}/${sample_name}
# $uhfdry_out stan.in \
#   && python3 finiteT.py -u 4 -g 31 --max 1 \
#   && python3 finiteT.py -u 8 -g 31 --max 3 \
#   && python3 finiteT.py -u 12 -g 31 --max 5 \
#   && python3 finiteT.py -u 24 -g 31 --max 7 \
#   :
# if [ $? -ne 0 ]; then
#   echo ${sample_name} failed.
#   failed_dirs="${failed_dirs} ${sample_name}"
#   nfails=`echo "$nfails + 1" | bc`
# else
#   echo ${sample_name} finished.
# fi

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
  echo ""
  echo "The following samples failed:"
  for failed_dir in $failed_dirs; do
    echo $failed_dir
  done
  exit 1
else
  echo ""
  echo "All samples finished."
  exit 0
fi

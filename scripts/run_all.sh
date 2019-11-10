#!/bin/bash

usage="$0 <folder> <method> <datatype> [reference]"

# method = algorithm8|jain_neal_split|triadic

folder=${1:? Usage: $usage}
method=${2:? "Usage: $usage"}
datatype=${3:? "Usage: $usage"}

if [ "$datatype" = "points3d" ]; then
  reference=${4:? "Usage: $usage"}
fi

# By default have the algorithm perform 10.000 MCMC steps
T=1000

mkdir -p logs

if [ "$datatype" = "points3d" ]; then

  for f in $folder/*.pnts.data; do
    fname="$(basename $f)"
    echo "(time  build/bin/noparama -d $f -c $datatype -T $T -a $method -r $reference) 2>&1 | tee logs/$fname.log"
    (time  build/bin/noparama -d $f -c $datatype -T $T -a $method -r $reference) 2>&1 | tee logs/$fname.log
  done

else
  
  for f in $folder/*.data.txt; do
    echo "build/bin/noparama -d $f -c $datatype -T $T -a $method"
    (time  build/bin/noparama -d $f -c $datatype -T $T -a $method) 2>&1 | tee logs/$f.log
  done

fi

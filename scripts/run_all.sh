#!/bin/sh

folder=${1:? "Usage: folder method, e.g. ../data/lines"}
method=${2:? "Usage: folder method, e.g. algorithm8|jain_neal_split|triadic"}

# By default have the algorithm perform 10.000 MCMC steps
T=10000

for f in $folder/*.data.txt; do
  echo "build/bin/noparama -d $f -c regression -T $T -a $method"
  time  build/bin/noparama -d $f -c regression -T $T -a $method
done

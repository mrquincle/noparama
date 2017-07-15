#!/bin/sh

folder=${1:? "Usage: folder method, e.g. folder=../data/lines"}
method=${2:? "Usage: folder method, e.g. method=algorithm8|jain_neal_split|triadic"}

T=10000

for f in $folder/*.data.txt; do
  echo "build/bin/noparama -d $f -c regression -T $T -a $method"
  build/bin/noparama -d $f -c regression -T $T -a $method
done

#!/bin/sh

directory=${1:? "directory"}

touch $directory/purity.txt
rm $directory/purity.txt
touch $directory/rand.txt
rm $directory/rand.txt
touch $directory/adjusted_rand.txt
rm $directory/adjusted_rand.txt

for f in $directory/*; do
  cat $f/results.score.txt | grep ^Purity | cut -f2 -d: >> $directory/purity.txt
  cat $f/results.score.txt | grep ^Rand | cut -f2 -d: >> $directory/rand.txt
  cat $f/results.score.txt | grep ^Adjusted | cut -f2 -d: >> $directory/adjusted_rand.txt
done

echo "You can find the results in purity.txt, rand.txt, and adjusted_rand.txt in $directory"
#cat $directory/purity.txt
#cat $directory/rand.txt
#cat $directory/adjusted_rand.txt

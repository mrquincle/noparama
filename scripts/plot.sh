#!/bin/sh

# if blah blah

cd ../build/bin

mv test_multivariate_normal_distribution.data ~/tmp

./test_multivariate_normal_distribution
gnuplot ../../plot/test_multivariate_normal_distribution.plot

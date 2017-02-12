#!/bin/bash

plot=${1:? "$0 test plot"}

cd ../build/bin

case $plot in 
  test_multivariate_normal_distribution|test_dirichlet)
    mv $plot.data ~/tmp
    ./$plot
    gnuplot ../../plot/$plot.plot
  ;;
esac


#!/bin/bash

plot_type=${1:? "$0 plot-type data-file"}
data_file=${2:? "$0 plot-type data-file"}

script_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
noparama_path="$script_path/.."
plot_path="$noparama_path/plot"

case $plot_type in 
test_multivariate_normal_distribution|test_dirichlet)
  gnuplot -e "filename=\"$data_file\"" "$plot_path/$plot_type.plot" 
  ;;
test_dirichlet_2d)
  head $data_file
  gnuplot -e "filename=\"$data_file\"" "$plot_path/test_multivariate_normal_distribution.plot" 
  ;;
*)
  echo "Unknown plot type"
  ;;
esac


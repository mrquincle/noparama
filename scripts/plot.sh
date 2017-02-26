#!/bin/bash

getopt --test > /dev/null
if [[ $? -ne 4 ]]; then
  echo "I’m sorry, `getopt --test` failed in this environment."
  exit 1
fi

SHORT=vp:d:
LONG=verbose,plot:,data:

PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [[ $? -ne 0 ]]; then
  exit 2
fi

eval set -- "$PARSED"

while true; do
    case "$1" in
      -v|--verbose)
	verbose=y
	shift
	;;
    -p|--plot)
	plot_type=$2
	shift 2
        ;;
    -d|--data) 
	data_file=$2
	echo $2
	shift 2
	;;
    --)
	shift
	break
	;;
     *)
	echo "Error in arguments."
	echo $2
	exit 3
	;;
    esac
done

usage() {
  echo "Usage: $0 --plot plot_type --data data_file"
}

if [ ! "$plot_type" ] || [ ! "$data_file" ]; then
  usage
  exit 4
fi

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


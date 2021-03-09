#!/bin/bash
# -----------------------------------------------------------------------------
# HLRS - Tensor optimization
#
# Author:          Johannes Gebert <gebert@hlrs.de>
# Created:         06.03.2021
# Last edit:       09.03.2021
# ------------------------------------------------------------------------------
# Check whether all requirements are met
# File checking done in bash. It's more natural doing it this way. Furthermore, it is not a part of the toolchain itself
# ------------------------------------------------------------------------------

stop=0

if [ ! -z $1 ]; then
	filename_in=$(basename -s ".csv" $1)
	filname_out_CR0=$filename_in"_CR0.csv"
	filname_out_CR1=$filename_in"_CR1.csv"
	filname_out_CR2=$filename_in"_CR2.csv"
	if [ ! -f $filename_out_CR0 ]; then
		stop=1
	fi
	if [ ! -f $filename_out_CR1 ]; then
		stop=1
	fi
	if [ ! -f $filename_out_CR2 ]; then
		stop=1
	fi
#	amnt_lines=$(cat $1 | wc -l)  Not clear how to handle it at the moment.
  if [ $stop == 0 ]; then
	    ./bin/tensor_optimizer_x86_64 $1 $amnt_lines
  else
      echo "Output file(s) already exist!"
  fi
else
	echo "Input file missing."
	echo "Usage: ./<this script> <input file>"
fi

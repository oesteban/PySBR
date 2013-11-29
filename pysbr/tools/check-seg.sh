#!/bin/sh 
#
# file       check_seg.sh
# brief      
# author     Gert Wollny (gert@die.upm.es)
# date       September, 2013
# ingroup    PySBR
#
# Evaluate the shape statistics for the input images 
# and save them
#
#
# Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban), and 
#                     gert@die.upm.es (Gert Wollny)
#                     with Biomedical Image Technology, UPM (BIT-UPM)
# All rights reserved.
# This file is part of PySBR.
#


if [ -f $HOME/.pysbr ]; then 
    root=$(cat "$HOME/.pysbr")
fi 
 

if [ "x$#" = "x1" ]; then 
    root="$1" 
fi 

if [ ! -d "$root" ] ; then 
    echo "the direcory '$root' isn't a directory, please specify a directory with the test data" 
    echo "on the command line of in the file '$HOME/.pysbr'"
    exit 1
fi 


sets="S040  S205  S289  S480  S549 S109  S214  S356  S487  S555 S116  S228"


maskflag="-K" 
maskfile="dat_phantom_seg_out.nii.gz"

if [ -f $root/${s}/bbox.nii.gz ]; then 
    echo "found mask file" 
    maskfile="bbox.nii.gz"
    maskflag="--inmask"
fi 

mkdir -p debug-results
for s in $sets; do
    for level in 0 1 2 3; do
        name="$root/${s}/dat_phantom_grad${level}_out.nii.gz"
        mask="$root/${s}/${maskfile}"
        echo $name
	if [ -f $name ]; then 
            ./sbr_get_area_stats.py -i $name $maskflag ${mask} -P 1.4 -B 0.7 -U 1.01  \
                --debug debug-results/${s}-${level}.nii.gz
	else 
	    echo "file '$name' not found"
	fi
    done
done 


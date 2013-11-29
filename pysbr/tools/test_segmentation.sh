#!/bin/sh 
#
# file       training_levels.sh
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

pushd $root 
sets=$(ls)
popd 

maskfile=bbox.nii.gz


bridge_ratio="0.7" 
mean_ratio="1.4"
max_volume="3200000" 
uphill="1.0"

outname=test

mkdir -p "$outname"


b=$bridge_ratio
m=$mean_ratio
v=$max_volume
u=$uphill

level=0 
                
for s in $sets; do
    echo $s
    if [ -d "$root/$s" ]; then 
        name="$root/${s}/dat_phantom_grad0_out.nii.gz"
        mask="$root/${s}/bbox.nii.gz"
        seg="$root/${s}/dat_phantom_seg_out.nii.gz"
    
        [ -e $outname/${s}_grade_${level}_spot.txt ] || \
            ./sbr_get_labelmask_and_thinned.py -i "$name" --inmask "${mask}" -U $u -P $m  -B $b -V $v  \
	    -m "$outname/${s}_grade_${level}_masks.nii.gz" -s "$outname/${s}_grade_${level}_spot.txt"
    
        ./sbr_dice.py -s "$outname/${s}_grade_${level}_masks.nii.gz" \
            -r "$root/${s}/dat_phantom_seg_out.nii.gz" >> "$outname/dice.txt" 
    fi
done
R < mean-singleclm.r --slave --args "$outname/dice.txt" 1 | sed -s "s/\[1\]//" | tee dice-stats.txt


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


sets="S063  S205  S289  S480  S549 S109  S214  S356  S487  S555 S116  S228"

maskfile=bbox.nii.gz


bridge_ratio="0.7 0.8 0.83 0.9 0.95 1.0" 
mean_ratio="0.0 0.7 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 2.0"
max_volume="3200000" 
uphill="1.0 1.1"

mkdir -p tune-path 

rm -f tune-path/dice-*.txt 

for b in $bridge_ratio; do 
    for m in $mean_ratio; do 
        for v in $max_volume; do 
            for u in $uphill; do 
                echo . 
                outname="${b}-${m}-${v}-${u}"
                mkdir -p "tune-path/$outname" 
                
                for s in $sets; do
                    
                    name="$root/${s}/dat_phantom_grad0_out.nii.gz"
                    mask="$root/${s}/bbox.nii.gz"
                    seg="$root/${s}/dat_phantom_seg_out.nii.gz"
                    
                    [ -e tune-path/$outname/${s}_grade_${level}_spot.txt ] || \
                        ./sbr_get_labelmask_and_thinned.py -i "$name" --inmask "${mask}" -U $u -P $m  -B $b -V $v  \
		            -m "tune-path/$outname/${s}_grade_${level}_masks.nii.gz" -s "tune-path/$outname/${s}_grade_${level}_spot.txt"
                    
                  ./sbr_dice.py -s "tune-path/$outname/${s}_grade_${level}_masks.nii.gz" -r "$root/${s}/dat_phantom_seg_out.nii.gz" >> "tune-path/dice-$outname.txt" 
                done
                value=$(R < mean-singleclm.r --slave --args "tune-path/dice-$outname.txt" 1 | sed -s "s/\[1\]//")
                echo "$outname $value" >> "tune-path/dice-mean.txt"

            done 
        done 
    done 
done
sort -nr -k 2 "tune-path/dice-mean.txt" >"tune-path/dice-mean-sorted.txt"

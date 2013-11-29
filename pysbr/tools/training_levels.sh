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


rm levels*.txt
rm dice_*.txt


maskfile=bbox.nii.gz

mkdir -p mask
for s in $sets; do
    for level in 0 1 2 3; do
        name="$root/${s}/dat_phantom_grad${level}_out.nii.gz"
        mask="$root/${s}/bbox.nii.gz"
        echo $name
	if [ -f $name ]; then 
            ./sbr_get_labelmask_and_thinned.py -i $name --inmask ${mask} -U 1.2 -P 0.0 -B 1.0  \
		-m mask/${s}_grade_${level}_masks.nii.gz -s mask/${s}_grade_${level}_spot.txt

            ./sbr_get_area_stats.py -i $name --inmask ${mask} -P 0.0 -B 1.0 -U 1.2 -a 400  >> levels_${level}.txt 
            
            if [ "$level" = "0" ]; then 
                ./sbr_dice.py -s mask/${s}_grade_${level}_masks.nii.gz -r "$root/${s}/dat_phantom_seg_out.nii.gz" >> dice_${level}.txt 
            fi 
	else 
	    echo "file '$name' not found"
	fi
    done
done 

for level in 0 1 2 3; do
    grep grav levels_${level}.txt | awk '{print $5 $6 $7 $8 $9 }' | \
        sed -e 's/[-a-z,]*=\([na0-9.]*\)/\1 /g' >levels${level}.txt
done 


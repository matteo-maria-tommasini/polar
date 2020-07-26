#!/bin/bash

MOL=$(find . -name '*.out.gz')

for zipped_gau in $MOL; do

   expanded_gau=$(basename ${zipped_gau/.out.gz/.out})
   dname=$(dirname $zipped_gau)
   name=${dname/.\//};
   echo $name

   gunzip --to-stdout $zipped_gau > $expanded_gau
 
   # spectra
   ../../polar.x -gau $expanded_gau -NP 3000 -min 0.0 -max 20.0 -g 0.35 > c84-$name.polar 
   ../../polar.x -r -gau $expanded_gau 
   ../../polar.x -d -gau $expanded_gau 
 
   # rotatory strengths
   cat $expanded_gau | sed -n '/Rotatory Strengths (R) in cgs/,/Rotatory Strengths (R) in cgs/p' | sed '1,2d' | sed '$d' | sed '$d' > R-$name.tmp 

   # states
   cat $expanded_gau | grep 'Excited State' | awk '{print $3,$4,$5,$6}' > states-$name.tmp

   paste states-$name.tmp R-$name.tmp > R-$name.dat
   rm -rf states-$name.tmp R-$name.tmp

   rm -rf $expanded_gau
done


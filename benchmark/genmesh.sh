#!/bin/bash
echo "Bash version ${BASH_VERSION}"

tableau=('0.05' '0.025' '0.0125' '0.00625')
#h = 0.05 0.025 0.0125
Nbssdom='32'

for i in ${tableau[*]}
do
    ./gmsh -2 -clmax $i -clmin $i -part $Nbssdom carre.geo -o carre_h${i}_N${Nbssdom}.msh
   # echo $i
done    

#echo $Nbssdom
echo "All done"


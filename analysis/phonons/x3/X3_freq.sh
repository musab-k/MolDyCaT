#!/bin/bash
cell_x=DUMMYA
cell_y=DUMMYA
cell_z=DUMMYA
#cd /Users/vhc08/Documents/GPP_Project/moldycat/
echo "#disp #E" > pho_X3_E_x.txt
for (( x=0; x<=10; x++ ))
do
  cd inputs
  disp=`echo "scale=2; -0.002+$x*0.0004" | bc`
  echo $disp
  ./gen_pos_X3.py $disp $cell_x $cell_y $cell_z
  sed "s/^CELL_DIMENSIONS.*$/CELL_DIMENSIONS $cell_x $cell_y $cell_z/" phonons_params.base | sed "s/^ATOMS_FILE_ADDRESS.*$/ATOMS_FILE_ADDRESS phonons_X3.txt/" > phonons_params.in
  cd ..
  ./moldycat phonons > temp
  E=`grep "System energy" temp | awk '{print $3}'`
  echo $disp $E >> pho_X3_E_x.txt
done

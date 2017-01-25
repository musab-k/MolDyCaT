#!/bin/bash
cell_x=3.55648
cell_y=3.55648
cell_z=3.55648
#cd /Users/vhc08/Documents/GPP_Project/moldycat/
echo "#disp #E" > pho_X1_E_x.txt
for (( x=0; x<=10; x++ ))
do
  cd inputs
  disp=`echo "scale=2; -0.002+$x*0.0004" | bc`
  echo $disp
  ./gen_pos_X1.py $disp
  sed "s/^CELL_DIMENSIONS.*$/CELL_DIMENSIONS $cell_x $cell_y $cell_z/" phonons_params.base | sed "s/^ATOMS_FILE_ADDRESS.*$/ATOMS_FILE_ADDRESS phonons_X1.txt/" > phonons_params.in
  cd ..
  ./moldycat phonons > temp
  E=`grep "System energy" temp | awk '{print $3}'`
  echo $disp $E >> pho_X1_E_x.txt
done
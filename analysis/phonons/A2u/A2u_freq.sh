#!/bin/bash
cell_x=CELLX
#cell_x=4.257
cell_y=CELLY
#cell_y=2.4577801
#cell_z=10
cell_z=3.35
#cd /Users/vhc08/Documents/GPP_Project/moldycat/
echo "#disp #E" > pho_A2u_E_x.txt
for (( x=0; x<=10; x++ ))
do
  cd inputs
  disp=`echo "scale=2; -0.004+$x*0.0008" | bc`
  echo $disp
  ./gen_pos_A2u.py $disp $cell_x $cell_y $cell_z
  sed "s/^CELL_DIMENSIONS.*$/CELL_DIMENSIONS $cell_x $cell_y $cell_z/" phonons_params.base | sed "s/^ATOMS_FILE_ADDRESS.*$/ATOMS_FILE_ADDRESS phonons_A2u.txt/" | sed "s/^MULTIPLY_CELL_DIRN.*$/MULTIPLY_CELL_DIRN 10 15 1/" > phonons_params.in
  cd ..
  ./moldycat phonons > temp
  E=`grep "System energy" temp | awk '{print $3}'`
  echo $disp $E >> pho_A2u_E_x.txt
done

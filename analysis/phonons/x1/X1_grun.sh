#!/bin/bash

gps='X1_omega.gnu'
grs='X1_freq.sh'
res='X1_gruResults.txt'

cell0=3.55648

x=-0.001
dx=0.001

echo "#x cell V omega" > $res 
for (( i=0; i<=2; i++ ))
    do
        cellnew=`echo $cell0 $x | awk '{print $1*(1+$2)}'`
        V=`echo $cellnew | awk '{print $1*$1*$1}'`
        
        sed "s:DUMMYA:$cellnew:" $grs > $grs.tmp2
        sed "s:DUMMYPHO:$x:" $grs.tmp2 > $grs.tmp

        sed "s:DUMMYA:$cellnew:" $gps > $gps.tmp2
        sed "s:DUMMYPHO:$x:" $gps.tmp2 > $gps.tmp

        chmod +x $grs.tmp
        ./$grs.tmp
        gnuplot $gps.tmp 2> $gps.ans
        omega=`tail -1 $gps.ans | awk '{print $1}'`

        echo $x $cellnew $V $omega >> $res

        x=`echo $x $dx | awk '{print $1+$2}'`
done

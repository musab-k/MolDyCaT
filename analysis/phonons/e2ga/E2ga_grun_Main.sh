#!/bin/bash

gps='E2ga_omegaFinder.gnu'
grs='E2ga_freq.sh'
res='gruResults_E2ga.txt'
phonon_results='pho_E2ga_E_x.txt'
grunScriptOrig='grunFinderGeneral.gnu'
grunScriptNew='E2ga_grunScript.gnu'

cell0x=4.257
cell0y=`echo $cell0x | awk '{print $1/sqrt(3)}'`
cell0z=3.35

x=-0.001
dx=0.001

echo "#x cell V omega" > $res 
for (( i=0; i<=2; i++ ))
    do
        cellnewx=`echo $cell0x $x | awk '{print $1*(1+$2)}'`
        cellnewy=`echo $cellnewx | awk '{print $1/sqrt(3)}'`
        V=`echo $cellnewx $cellnewy $cell0z | awk '{print $1*$2*$3}'`
        
        sed "s:CELLX:$cellnewx:" $grs > $grs.tmp2
        sed "s:CELLY:$cellnewy:" $grs.tmp2 > $grs.tmp3
        sed "s:DUMMYPHO:$phonon_results:" $grs.tmp3 > $grs.tmp

        sed "s:DUMMYA:$cellnewx:" $gps > $gps.tmp2
        sed "s:DUMMYPHO:$phonon_results:" $gps.tmp2 > $gps.tmp

        chmod +x $grs.tmp
        ./$grs.tmp
        gnuplot $gps.tmp 2> $gps.ans
        omega=`tail -1 $gps.ans | awk '{print $1}'`

        echo $x $cellnewx $V $omega >> $res

        if [ "$x" = "0" ]; then
            omega_zero=`echo $omega | awk '{print $1}'`
        fi

        x=`echo $x $dx | awk '{print $1+$2}'`
done

sed "s:FILE:$res:" $grunScriptOrig > $grunScriptNew.tmp2
sed "s:OMEGA:$omega_zero:" $grunScriptNew.tmp2 > $grunScriptNew.tmp
gnuplot $grunScriptNew.tmp

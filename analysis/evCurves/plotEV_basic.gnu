#example gnuplot script

#required for angstroms
set encoding iso_8859_1
###Output to pdf. Will use different colours for lines:
set term pdf enhanced font "Times New Roman, 9" 
set output "evCurves_basic.pdf"

###Output to eps. These can easily be included in a TeX document. They will be black and white. 
#set term postscript eps enhanced 
#set output "Example.eps"

set size ratio 1.2

set nokey
#set key spacing 1.5
set autoscale
set xrange[1.0:1.8]
set yrange[-9:-4]
unset log
unset label
set xtic auto
set xtics 0.2
set ytic auto

#set title "All E-V curves"
set xlabel "Nearest neighbour separation (\305)"
set ylabel "Energy (eV/atom)"
set label "Diamond" at 1.6,-8.5
set label "Graphite" at 1.32,-8.6
set label "Linear Chain" at 1.41,-6.6

plot 'energy_vol_Graphene.txt' w l title "graphite", 'energy_vol_diamond512.txt' w l title "diamond",'energy_vol_linear_chain.txt' w l title "linear chain"
 # removed nanotubes: 'energy_vol_nanotube.txt' w l title "nanotube", 

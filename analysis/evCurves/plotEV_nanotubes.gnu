#example gnuplot script

#required for angstroms
set encoding iso_8859_1
###Output to pdf. Will use different colours for lines:
set term pdf enhanced font "Times-Roman, 9" 
set output "evCurves_nanotubes.pdf"

###Output to eps. These can easily be included in a TeX document. They will be black and white. 
#set term postscript eps enhanced 
#set output "Example.eps"

#set nokey
set key top right
#set key spacing 1.5
set autoscale
#set xrange[1.23:1.76]
#set yrange[-10:0]
unset log
unset label
set xtic auto
set ytic auto

#set title "All E-V curves"
set xlabel "Nearest neighbour separation (\305)"
set ylabel "Energy (eV/atom)"

plot 'energy_vol_nano_3_0.txt' w l title "(3,0)",  'energy_vol_nano_3_3.txt' w l title "(3,3)",  'energy_vol_nano_13_0.txt' w l title "(13,0)"

#example gnuplot script

#required for angstroms
set encoding iso_8859_1
###Output to pdf. Will use different colours for lines:
set term pdf enhanced font "Times-Roman, 9" 
set output "evCurves_conv_armchair_nanotube.pdf"

###Output to eps. These can easily be included in a TeX document. They will be black and white. 
#set term postscript eps enhanced 
#set output "Example.eps"

set nokey
#set key spacing 1.5
set autoscale
#set xrange[1.2:1.8]
#set yrange[-9:-4]
unset log
unset label
set xtic auto
set ytic auto

#set title "All E-V curves"
set xlabel "Number of Atoms"
set ylabel "Energy (eV/atom)"

plot 'energy_conv_nanotube.txt' w l title "Armchair carbon nanotube"

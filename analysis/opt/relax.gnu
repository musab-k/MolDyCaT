#example gnuplot script

#required for angstroms
set encoding iso_8859_1

###Output to pdf. Will use different colours for lines:
set term pdf enhanced font "Times-Roman, 6"
set output "relaxers.pdf"

###Output to eps. These can easily be included in a TeX document. They will be black and white. 
#set term postscript eps enhanced 
#set output "Example.eps"

set key top right
set key spacing 1.5
set autoscale
#set xrange[0:50]
#set yrange[0:1]
unset log
unset label
set xtic auto
set ytic auto

set xlabel "Relaxation step"
set ylabel "System energynergy (eV/\305^3)"
plot "sd.out" u 1:3 ti "SD", "cg.out" u 1:3 ti "CG", "dmdum2.out" u 1:3 ti "DMD-UM", "dmdcm.out" u 1:3 ti "DMD-CM"

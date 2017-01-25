#example gnuplot script

#required for angstroms
set encoding iso_8859_1

###Output to pdf. Will use different colours for lines:
set term pdf enhanced font "Times-Roman, 9"
set output "grapheneC11C21.pdf"

###Output to eps. These can easily be included in a TeX document. They will be black and white. 
#set term postscript eps enhanced 
#set output "Example.eps"

set nokey
#set key top right
#set key spacing 1.5
#set autoscale
#set xrange[0:50]
#set yrange[0:1]
unset log
unset label
set xtic auto
set ytic auto
set pointsize 2

set xlabel "x"
set ylabel "Energy per unit volume (eV/\305^3)"

a=1;b=1;c=1;
f(x)=a*x**2+b*x+c
fit f(x) 'grapheneC11C21_240.txt' u 1:2 via a,b,c
plot 'grapheneC11C21_240.txt' u 1:2 lt -1 t 'Data', f(x) lt 3 t 'Fit'
print "" 
print "" 
print "" 
print "C11-C21 = " 
print 1.602177*a

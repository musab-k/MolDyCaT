#example gnuplot script

#required for angstroms
set encoding iso_8859_1

###Output to pdf. Will use different colours for lines:
set term pdf enhanced font "Times-Roman, 6"
set output "diamondC44.pdf"

###Output to eps. These can easily be included in a TeX document. They will be black and white. 
#set term postscript eps enhanced 
#set output "diamondC44.eps"

set key top right
set key spacing 1.5
set autoscale
#set xrange[0:50]
#set yrange[0:1]
unset log
unset label
set xtic auto
set ytic auto

set xlabel "x"
set ylabel "Energy per unit volume (eV/\305^3)"

a=1;b=1;c=1;                                                                                                       
f(x)=0.5*a*x**2+b*x+c
a=1;b=1;c=1;
fit f(x) 'diamondC44_tama_432.txt' via a,b,c
plot 'diamondC44_tama_432.txt' t 'Data', f(x) t 'Fit'

print ""
print ""
print ""
print "C44 = " 
print 1.602177*a

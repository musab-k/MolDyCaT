a=1;b=1;c=1;
f(x)=0.5*a*x**2+b*x+c
a=1;b=1;c=1;
fit f(x) 'DUMMY' via a,b,c
print 1.602177*a
plot 'DUMMY', f(x)

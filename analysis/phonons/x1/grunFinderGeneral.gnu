f(x)=a*x+b
a=1; b=1;
fit f(x) 'FILE' u 3:4 via a,b
V=44.9843
w=OMEGA
g=-a*V/w
print g

plot 'FILE' u 3:4, f(x)

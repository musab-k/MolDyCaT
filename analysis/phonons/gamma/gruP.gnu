f(x)=a*x+b
a=1; b=1;
fit f(x) 'gruResults.txt' u 3:4 via a,b
V=44.9843
w=36.7897170218831
g=-a*V/w
print g

plot 'gruResults.txt' u 3:4, f(x)

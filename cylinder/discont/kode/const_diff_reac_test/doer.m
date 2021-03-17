a = csvread("numsol.csv");
xize = size(a)/100
xize = xize(1)
_z = a(1:100,1:100);
m = min(min(_z))
M = max(max(_z))
n = 100;

for i=0:(xize-1)
j=100*i
z = a(j+1:j+100,1:100);
#z = a(201:300,1:100)*100;
z = z - m;
z = z / M;
z = n * z;
figure(i+1);
clf;
##image(z);
##colormap (spring (n));
##filename = sprintf("fig%d.png",i)
##print(i+1,filename);
contour(z, linspace(0, n, 25));
endfor
a = csvread("numsol.csv");
xize = size(a)/100
xize = xize(1)
z = a(1:100,1:100);
m = min(min(z))
M = max(max(z))
n = 100;
cmap = colormap (summer (100));

for i=0:(xize-1)
j=100*i
z = a(j+1:j+100,1:100);
#z = a(201:300,1:100)*100;
#z = z - m;
#z = z / M;
#z = n * z;
figure(i+1);
##clf;
##image(z);
##colormap (summer (100));
##colorbar
##filename = sprintf("col%d.png",i)
##print(i+1,filename);
clf;
contourf(z);
caxis([m M]);
colormap(cmap);
colorbar
xlabel("R", "fontsize", 15);
ylabel("Z", "fontsize", 15);
filename = sprintf("con%d.png",i);
print(i+1,filename);
endfor

a = csvread("uniformity.csv");
image(a);
#caxis([m M]);
colormap(cmap);
colorbar
print -dpng uniformity.png
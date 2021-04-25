a = csvread("chnumsol.csv");
##aa = csvread("zconst_diff_ex.csv");
xize = size(a)/100
xize = xize(1)
z = a(1:100,1:100);
x = 1:100;
m = min([0 min(min(z))])
M = max(max(z))
n = 100;
#cmap = colormap (summer (100));
##cmap = colormap (summer);

for i=0:10#(xize-1)
j=100*i
z = a(j+1:j+100, :);
##zz = aa(i+1, :);
#z = a(201:300,1:100)*100;
#z = z - m;
#z = z / M;
#z = n * z;
figure(i+1);
##colormap(cmap);
##caxis([0 1]);
##plot(x, z(:,1)
###  , "linestyle", "none", "marker", "x"
##  , x, zz
###  , "linestyle", "none", "marker", "x"
## );
###ylim([-0.5, 1]);
##grid on;
##legend(["numerical";"exact"] ,"fontsize",15);
##clf;
##image(z);
##colorbar
##filename = sprintf("col%d.png",i)
##print(i+1,filename);
clf;
contourf(z);
##caxis([0 2]);
##image(z);
colorbar
xlabel("R", "fontsize", 15);
ylabel("Z", "fontsize", 15);
title(num2str(i));
##saveas(sprintf("ch%d.png",i));
filename = sprintf("chuni%d.png",i);
print(i+1,filename);
endfor

##a = csvread("uniformity.csv");
###image(a);
##contourf(a);
###caxis([m M]);
##caxis([2 3]);
##colormap(cmap);
##colorbar
###print -dpng uniformity.png








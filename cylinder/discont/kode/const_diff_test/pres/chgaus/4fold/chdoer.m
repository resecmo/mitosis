a = csvread("chnumsol.csv");
##aa = csvread("zconst_diff_ex.csv");
xize = size(a)/50
xize = xize(1)
z = a(1:50,1:end);
x = 1:50;
m = min([0 min(min(z))])
M = max(max(z))
n = 50;
cmap = colormap (viridis (1000));
##cmap = colormap (summer);
phi = linspace(-1,1,100)*pi/2;
circr = 5*cos(phi);
circz = 26 + 5*sin(phi);
for i=0:(xize-1)
j=50*i
z = a(j+1:j+50, :);
##zz = aa(i+1, :);
#z = a(201:300,1:100)*100;
#z = z - m;
#z = z / M;
#z = n * z;
figure(i+1);
##colormap(cmap);
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
colormap(cmap);
##caxis([0, max(max(z))+0.001]);
hold on; plot(circr,circz, 'color', 'red'); hold off;
##caxis([0 2]);
##image(z);
colorbar
xlabel("R", "fontsize", 15);
ylabel("Z", "fontsize", 15);
title(num2str(i*200));
##saveas(sprintf("ch%d.png",i));
filename = sprintf("gus%02d.png",i);
##print(i+1,filename);
endfor

##a = csvread("uniformity.csv");
###image(a);
##contourf(a);
###caxis([m M]);
##caxis([2 3]);
##colormap(cmap);
##colorbar
###print -dpng uniformity.png








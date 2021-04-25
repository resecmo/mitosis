a = csvread("rnumsol.csv");
aa = csvread("rconst_diff_ex.csv");
xize = size(a)/100
xize = xize(1)
z = a(1:100,1:100);
x = (1:100)/100;
m = min(min(z))
M = max(max(z))
n = 100;
cmap = colormap (summer (100));

for i=0:(xize-1)
j=100*i
z = a(j+1:j+100,1:100);
zz = aa(i+1, 1:100);
abser = abs(z-zz);
err = min(abser, abser./abs(zz));
#z = a(201:300,1:100)*100;
#z = z - m;
#z = z / M;
#z = n * z;
figure(i+1);
clf;
colormap(cmap);
caxis([m M]);
subplot(1,2,1);
  plot(x, z(10,:), x, zz);
  ylim([-0.5, 1]);
  grid on;
  xlabel('r', 'fontsize', 15);
  ylabel('[a*]', 'fontsize', 15);
  title(sprintf('%d^{th} timestep', i), 'fontsize', 15);
  legend(['numerical';'exact'], 'fontsize', 15);
subplot(1,2,2);
  plot(x,err(10,:) * 100);
  ylim([0, 0.2] * 100);
  grid on;
  xlabel('r', 'fontsize', 15);
  ylabel('min(\delta, \epsilon)', 'fontsize', 15);
  title('err, %', 'fontsize', 15);
saveas(sprintf("r%d.png", i));
##clf;
##image(z);
##colorbar
##filename = sprintf("col%d.png",i)
##print(i+1,filename);
##clf;
##contourf(z);
##colorbar
##xlabel("R", "fontsize", 15);
##ylabel("Z", "fontsize", 15);
##filename = sprintf("con%d.png",i);
##print(i+1,filename);
endfor

##a = csvread("uniformity.csv");
###image(a);
##contourf(a);
###caxis([m M]);
##caxis([2 3]);
##colormap(cmap);
##colorbar
#print -dpng uniformity.png








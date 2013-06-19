D = dlmread('pdegdist.txt');
x = D(:,1);
y = D(:,2);
for i = 1:30
	xax(i) = x(i);
	yax(i) = y(i);
end

scatter(xax,yax);
t = 0:1:30;
hold on;
plot(t, 1/4*(3/4).^(t-3),'r');

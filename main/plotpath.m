P = dlmread('ppath.txt');
x = P(:,1);
y = P(:,2);
logy = log(y);

scatter(x,y);

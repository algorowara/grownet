M = dlmread('pgrowdata.txt');
x = M(:,1);
y = M(:,2);
age = M(:,3);

scatter3(x, y, age, '*');

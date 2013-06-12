M = dlmread('points.txt');
x = M(:, 1);
y = M(:, 2);
z = M(:, 3);
scatter3(x, y, z);

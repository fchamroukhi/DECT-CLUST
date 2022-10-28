function [X1,X2] = grid2d(m1,m2)
X1_axis=linspace(0,1,m1)';
X2_axis=linspace(0,1,m2)';

[X1,X2] = meshgrid(X1_axis,X2_axis);
X1 = X1(:);  
X2 = X2(:);

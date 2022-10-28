function solution = Bsplinefit(y,x,knots,M)
%
% [yfit_spline, coeff_spline, deg_freedom_spline, reg_matrix_spline, hatmatrix_spline] = Bspline(y,x,knots,M)
% fit a Bspline of order M to a noisy data y
% Inputs: 
%
%       y: the data [nx1] column vector
%       x: predictors [nx1] column vector
%       M: spline order (polynomial degree = M-1)
%       knots: interior knots [Kx1] ([1xK]) (K the number of interior knots)
%
% Outputs: 
%       
%       yfit: the fitted function [nx1]
%       coeff: the splines coefficients [(K+M)x1] column vector
%       deg_freedom: degree of freedom of the model = trace(hatmatrix)
%       reg_matrix: Xij = [h_1(xi),..,h_j(xi),...,h_{M+K}(xi)]regression matrix of the model [nx(K+M)]
%       hatmatrix: H=X(X'X)^{-1}X' the hat matrix
% (yfit = X*coeff = hatmatrix*y)
%
% Faicel Chamroukhi October 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(x,2)~=1,   x=x'; end
if size(y,2)~=1,   y=y'; end
if size(knots,1)~=1,   knots=knots'; end
if length(x)~=length(y),   error('x and y must have the same length');end

n = length(y);

knots = [x(1) knots x(end)+x(end)-x(end-1)];% interior and the two boudaries knots
% knots = [x(1) knots x(end)];% interior and the two boudaries knots
X = bsplinebasis(x,knots,M);
% plot(X,'x-')
% pause
 coeff = inv(X'*X)*X'*y;
%coeff =  (X'*X)\(X'*y);

yfit =  X*coeff;

sigma2 = (1/n)*sum((y-yfit).^2); %variance

reg_matrix = X;
% hatmatrix = X*inv(X'*X)*X';
hatmatrix = X*((X'*X)\X');

deg_freedom = trace(hatmatrix);

solution.yfit = yfit;
solution.coeff = coeff;
solution.sigma2 = sigma2;
solution.deg_freedom = deg_freedom;
solution.reg_matrix = reg_matrix;
solution.hatmatrix = hatmatrix;






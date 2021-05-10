function X = naturalsplinebasis(x,knots,M)
%
% X = naturalsplinebasis(x,knots,M)
% construct the regression matrix  for a natural cubic spline regression model of order M
%
% Inputs: 
%
%       x: predictors [nx1] column vector
%       M: spline order (polynomial degree = M-1)
%       knots: interior knots [Kx1] ([1xK]) (K the number of interior knots)
%
% Outputs: 
%
%       X: regression matrix of the model such that: Xij = [h_1(xi),..,h_j(xi),...,h_{M+K}(xi)] [nx(K+M)]
%
% Faicel Chamroukhi October 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = M-1;
if size(x,2)~=1,    x=x'; end
X=[];
% for ord = 0:1
%     X =[X x.^(ord)];%[1 t t.^2 t.^3 t.^p;......;...]
% end
X =[ones(length(x),1) x];
% X = splinebasis(x,knots,M);
% X(:,[2 3])=[];
% knots = knots;
knots = [x(1) knots x(end)];%!!!!!!!!! NB
K = length(knots);
xi_K_1 = knots(K-1);
xi_K = knots(K);

dK_1 = ((max(x-xi_K_1,0)).^3-(max(x-xi_K,0)).^3)/(xi_K - xi_K_1);

for k=1:K-2 
    xi_k = knots(k);
    dk = ((max(x-xi_k,0)).^p-(max(x-xi_K,0)).^p)/(xi_K - xi_k);
    X =[X dk-dK_1];
end
function X = splinebasis(x,knots,M)
% X = splinebasis(x,knots,M)
% construct the spline regression matrix  for spline regression model of
% order M
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
p=M-1; % polynomial degree
if size(x,2)~=1, x=x'; end

X=[];
for ord = 0:p
    X =[X x.^(ord)];%[1 t t.^2 t.^3 t.^p;......;...]
end
%
K = length(knots);
for k=1:K
    X =[X max(x-knots(k),0).^p];
end


% p=M-1; % polynomial degree
% if size(x,2)~=1, x=x'; end
% 
% X=[];
% for ord = 0:p
%     X =[X x.^(ord)];%[1 t t.^2 t.^3 t.^p;......;...]
% end
% 
% 
% close, subplot(1,2,1), plot(X,'x-')
% %
% K = length(knots);
% Y = [];
% for k=1:K
%     Y =[Y max(x-knots(k),0).^p];
% end
% 
% subplot(1,2,2), plot(X,'x-')
% pause
function B = bsplinebasis(x, t, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% B = bsplinebasis(x, t, M) : construct the spline regression matrix  for Bspline regression model of
% order M
%
% Inputs:
%
%      x: Vector of points at which to evaluate the b-spline basis
%      t: Vector of knots.  Normally the first and last will be repeated m times
%         here t represents the knots sequence (including the two boundary knots)
%      M: order of B-spline (degree of poynomial pieces p = m-1)
% Outputs: 
%  
%      B: the B-spline Basis (yfit = B*coeff)
%
%
%
%
%
%
%
% % t is tau in the document
%Faicel Chamroukhi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if size(t,1)~=1 , t=t';end % make a row vector

m = length(x);
% repeat the first and the last knots m times
t = [t(1)*ones(1,M-1) t t(end)*ones(1,M-1)];
%j th basis function : 1 <= j <= length(t) - M -1 ; There are K + 2*M knots ; M = length(t) - 1;
B = zeros(m, length(t) - M);

for j=1:length(t) - M %is the same as j=1 : K+M, K is the number of interior knots
    B(:,j)=bsplineBasis_j(x, t, j, M);% j-1 0 <= j <= length(t) - n - 2.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Bj = bsplineBasis_j(x, t, j, M)
% bsplineBasis: compute a bspline basis function of order m.
%
%    Bj = bsplineBasis(x, t, j, m)
%
% Evaluates the j^th B-spline basis function of degree m for a given set of
% knots t at points x using the Cox-de Boor recursion formula.
%
% Inputs:
%
%   x: Vector of points at which to evaluate the b-spline basis.
%   t: Vector of knots.  Normally the first and last will be repeated n times.
%   j: Scalar specifying which basis function from the basis function
%   set.  1 <= j <= length(t) - m - 1.
%   m: Order of basis functions.  Default m = 4 (cubics since degree= m-1 = 3).
%
% Outputs:
%
%   Bj: Basis function evaluated at the points in x.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Default parameter setting.
  if(nargin < 4),  M = 4;  p=3;  end% by defalut degree = 3 (cubic splines)  
  p = M-1;% polynomial degree
  
  % There are K knots. % j=1:L = K+M
  L = length(t); % L=K+M

  % Check validity of j.
  if((j < 1) || (j > L - M))
    error([ 'Parameter j = ' num2str(j) ' not in range [ 1, ' num2str(L - M) ' ]']);
  end
  
  %% construct the b-spline recursively.
  if(p == 0) % piecewise constatnt
    % Base case.
    %if(j+2 == length(t))
    %  y = ((x >= t(j+1)) & (x <= t(j+2)));
    %else
      Bj = ((x >= t(j)) & (x < t(j+1)));
    %end
   
  else
    % If the two knots in the denominator are equal, we ignore the term (to not devide by zero).
    denom1 = t(j + M - 1) - t(j);
    if(denom1 == 0)
      Bj = zeros(size(x));
    else
      Bj = (x - t(j)) .* bsplineBasis_j(x, t, j, M-1) / denom1;
    end
    
    denom2 = t(j + M) - t(j + 1);
    if(denom2 ~= 0)
      Bj = Bj + (t(j + M) - x) .* bsplineBasis_j(x, t, j+1, M-1) / denom2;
    end
  end
  

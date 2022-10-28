function [NV, B] = designSRMF(V, X, mixingOption, regressionOptions)
% Construction of the regression desing matrices
% FC
% NV : matrices of neighberhoods in the coordinates matrix V

neighb = 1;
    
n = length(V);
%NV = cell{n,1};
for i=1:n
    si1 = V(i,1); si2 = V(i,2); % in 2d
    Sj1 = V(:,1); Sj2 = V(:,2);
    delta_ij = max(abs(si1 - Sj1),abs(si2 - Sj2));
    
    Si = find(delta_ij  <= neighb);
    %Si (Si == i) = [];
    
    V_i = V(Si,:);
    NV{i} = V_i;
%         size(V_i)
%         pause

end


% construct the design matrix for a polynomial regression of degree p or a spline or B-spline of order M
switch regressionOptions.basis
    case 'polynomial'
        p = regressionOptions.p;
        B = designmatrix_Poly_Reg(X,p);
    case 'spline'
        % uniform knots locations %
        knots = linspace(X(1),X(end),regressionOptions.nknots + 2); % including the first and the two boudaries knots
        % knots(end-1) = x(end)+x(end)-x(end-1)+x(end)-x(end-1);
        %knots(end) = X(end)+X(end)-X(end-1);
        % knots = [x(1) knots(2:end-1) x(end)+x(end)-x(end-1)];% interior and the two boudaries knots
        M = regressionOptions.spline_order;
        dimBeta =  M + regressionOptions.nknots;
        %
        knots =  knots(2:end-1);    % take the interior knots
        B = splinebasis(X, knots, regressionOptions.spline_order);
        
    case'B-spline'
        % uniform knots locations %
        knots = linspace(X(1),X(end),regressionOptions.nknots+2); % including the first and the two boudaries knots
        % knots(end-1) = x(end)+x(end)-x(end-1)+x(end)-x(end-1);
        knots(end) = X(end)+X(end)-X(end-1);
        % knots = [x(1) knots(2:end-1) x(end)+x(end)-x(end-1)];% interior and the two boudaries knots
        M = regressionOptions.Bspline_order;
        dimBeta =  M + regressionOptions.nknots;
        
        %knots =  knots(2:end-1);    % take the interior knots
        
        B = bsplinebasis(X, knots, regressionOptions.Bspline_order);
    otherwise
        error('not included regression model');
end

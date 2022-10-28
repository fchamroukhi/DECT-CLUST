function [V, B] = designSRM(V, X, mixingOption, regressionOptions)
% Construction of the regression desing matrices
% FC

if strcmp(mixingOption, 'softmax')
    V = [ones(size(V,1), 1) V ];
    %V = [ones(size(V,1), 1) V V(:,1).^2 V(:,1).*V(:,2) V(:,2).^2];
    %V = [ones(size(V,1), 1) V V(:,1).^2 V(:,1).*V(:,2) V(:,2).^2 V(:,1).^3 (V(:,1).^2).*V(:,2) V(:,1).*(V(:,2).^2) V(:,1).^3];
    %V = [ones(size(V,1), 1) V V(:,1).^2 V(:,1).*V(:,2) V(:,2).^2 V(:,1).^3 (V(:,1).^2).*V(:,2) V(:,1).*(V(:,2).^2) V(:,1).^3 V(:,1).^4 (V(:,1).^3).*V(:,2)    (V(:,1).^2).*(V(:,2).^2)    V(:,1).*(V(:,2).^3)  V(:,2).^4];
    
    %  [n, m] = size(Y);
    %  neig = 2;
    %  Y_tmp = [zeros(neig,m); Y ; zeros(neig,m)];
    %
    %  V = zeros(n, m+1+2+3);
    %  for i=neig+1:n
    %      mi = mean(Y_tmp(i-neig: i +neig, :));
    % V(i,:) = [1 Vin(i,:) V(i,1)^2 V(i,1)*V(i,2) V(i,2)^2 mi];
    %  end
    % % Bstack = [];
    % %  for i=1:n
    % %      Bstack = [Bstack; [ones(m, 1)*V(i,:) B]];
    % %  end
    % %  dimBeta = size(Bstack, 2);
    % %  B = Bstack(1:m,:);
else %Gaussian
    V;
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

function [V, B] = designSRM_mrf(V, X, Y, mixingOption, regressionOptions)
% Construction of the regression desing matrices
% FC

if strcmp(mixingOption, 'softmax')
 [n, m]= size(Y);
 coord = V;
    V =  V/max(max(V));
    %V = [ones(size(V,1), 1) V ];
    %V = [ones(size(V,1), 1) V V(:,1).^2 V(:,1).*V(:,2) V(:,2).^2];
    %V = [ones(size(V,1), 1) V V(:,1).^2 V(:,1).*V(:,2) V(:,2).^2 V(:,1).^3 (V(:,1).^2).*V(:,2) V(:,1).*(V(:,2).^2) V(:,1).^3];
    %V = [ones(size(V,1), 1) V V(:,1).^2 V(:,1).*V(:,2) V(:,2).^2 V(:,1).^3 (V(:,1).^2).*V(:,2) V(:,1).*(V(:,2).^2) V(:,1).^3 V(:,1).^4 (V(:,1).^3).*V(:,2)    (V(:,1).^2).*(V(:,2).^2)    V(:,1).*(V(:,2).^3)  V(:,2).^4];
  
    
    
    obs = Y;
    
    neigb = 100;

C = zeros(n, m);
for i=1:n
    si1 = coord(i,1); si2 = coord(i,2);
    delta_ij = zeros(n,1);
    for j=1:n
        sj1 = coord(j,1); sj2 = coord(j,2);
        delta_ij(j) = max(abs(si1 - sj1),abs(si2 - sj2));
    end
    
    Si = find(delta_ij  <= neigb);
    Si (Si == i) = [];
    
    cardSi = length(Si);
    

    Cs = obs(Si, :);
    
    %[sum(Cs == obs(i)) length(Cs)]
    
    %muCs_k = zers(1, K);
    %for k=1:K
    %    muCs_k(k) = sum(Cs == k)/cardSi;
    %end
    
    muCs_i = mean(Cs); %sum(Cs == obs(i))/cardSi;
    
    %delta_ij(i,i) = [];
    
    C(i, :) = muCs_i;
end
V = [ones(n,1) C];
  
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
% neighb = 2;  
%     n = length(V);
% C = zeros(n, 4);
% for i=1:n
%     si1 = V(i,1); si2 = V(i,2); % in 2d
%     %     delta_ij = zeros(n,1);
%     %         for j=1:n
%     %             sj1 = coord(j,1); sj2 = coord(j,2);
%     %             delta_ij(j) = max(abs(si1 - sj1),abs(si2 - sj2));
%     %         end
%     Sj1 = V(:,1); Sj2 = V(:,2);
%     delta_ij = max(abs(si1 - Sj1),abs(si2 - Sj2));
%     
%     Si = find(delta_ij  <= neighb);
%     %Si (Si == i) = [];
%     
%     muCs_i = mean(V(Si,:));
%     C(i, :) = [V(i,:) muCs_i];
% end
% V = C;   
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

function param = Initialize_SRMF(V, B, Y, K, mixingOption, init_kmeans, try_algo)
%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%%%%%%%%%%%%%%%%%% FC
%V = data.VoxelCoordinates;
%% Energy Curves
%X = data.WavelengthLevels;
%Y = data. ReflectanceValues;

[n, m]=size(Y);

q1 = size(V,2);
%basis = regressionOption.basis;
%X = regressionOption.design;

%%   1. Initialize the mixing weights function's parameters
if strcmp(mixingOption, 'softmax')
    if try_algo ==1; param.Alphak = zeros(q1,K-1);else, param.Alphak = rand(q1,K-1);end    
else
    %% Gaussian gating
    
    % Gating Net parameters
    sol_km = myKmeans(Y,K);
    klas = sol_km.klas;
    %param.Mus = sol_km.muk;
    
    [n, p] = size(V);
    param.Mus = zeros(K, p);
    param.R = zeros(p, p, K);
    param.Theta = zeros(p, p, K);
    param.Alphak = zeros(1,K);
   
    rndind = randperm(n);
    for k=1:K
        Vk = V(klas==k, :);
        nk = size(Vk,1);
        param.Alphak(k) = nk/n;

        %param.Mus(k,:) = V(rndind(k),:) ;%
        param.Mus(k,:) = mean(Vk);
        z = (Vk - ones(nk,1)*param.Mus(k,:));
         param.R(:,:,k)  = z'*z/nk;%eye(p);%cov(Xk);
        %param.R(:,:,k)  = eye(p);%cov(Xk);
        
%sk = sum(sum(z.*z))/(m*nk);
%param.R(:,:,k) = sk*eye(p);
        
        
        %z'*z/nk;
        %     if p<n
        %         Theta(:,:,k) = inv(R(:,:,k));
        %     else
        %param.Theta(:,:,k) = eye(p);
        %     end
        
        %param.R(:,:,k)  = eye(p)*sum(sum(z.*z))/(m*nk);
    end
end

%%2. Initialize the regression parameters betak and the variances sigma2k
dimBeta = size(B, 2);

Betak   = zeros(dimBeta,K);
Sigmak2  = zeros(K,1);
%Sigmak2  = zeros(m, m, K);

% 2. betak and sigma2k
if init_kmeans
    max_iter_kmeans = 400;
    n_tries_kmeans = 5;
    verbose_kmeans = 0;
    
    sol_kmeans = myKmeans(Y, K,n_tries_kmeans, max_iter_kmeans, verbose_kmeans);
    klas = sol_kmeans.klas;
else % random partition
    klas = zeros(n,1);
    for i=1:n; klas(i) = randperm(K,1); end
end

% Ytild = reshape(Y',[],1);
% Btild = repmat(B,n,1);
% 
% Zik = (klas*ones(1,K))==(ones(n,1)*[1:K]);
for k=1:K
    Yk = Y(klas==k ,:); %if kmeans
    [nk, m] = size(Yk);
    
    Btild = repmat(B,nk,1);
    
    
    %[param_init reg_matrix] =  init_regression_param(Yk, t, spline_type, spline_order, knots, try_algo);
    Yktild = reshape(Yk',[],1);
    % coeff = inv(Xtild'*Xtild)*Xtild'*Ytild;
    
    betak =  (Btild'*Btild)\(Btild'*Yktild);
    Betak(:,k) = betak;
    
    Ykfit =  B*betak;
    
    sigmak2 = sum(sum((Yk-ones(nk,1)*Ykfit').^2))/(m*nk);
    
    
    Sigmak2(k)= sigmak2;
    %Sigmak2(:,:,k) = cov(Yk);
    %             Sigmak2(:,:,k) = diag(diag(Sigmak2(:,:,k)));
%%
%             zik = Zik(:,k);
%             % update of the regression coefficients
%             tmp =  repmat(zik,1,m);% [m x n]
%             Wk = reshape(tmp',[],1);%cluster_weights(:)% [mn x 1]
%             % same thing
%             % temp =  repmat(tauik,1,m)';% [m x n]
%             % cluster_weights = cluster_weights(:);
%             wYk = Wk.*Ytild; %cluster k
%             wBk = (Wk*ones(1,dimBeta)).*Btild;%[(n*m)*(M+nknots)]
%             % maximization w.r.t betak: Weighted least squares
%             %betak  = inv(phik'*phik + 1e-4*eye(dimBeta))*phik'*Yk;
%             betak  = (wBk'*wBk)\(wBk'*wYk);
%             Betak(:,k) = betak;
%             %B*betak;
%             % update the variance
%             sigmak2 = sum((wYk - wBk*betak).^2)/sum(Wk);%(m*sum(tauik));%
%             %             % sigmak2 = (1-gama)*sigmak2 + gama*Q;
%                          Sigmak2(k) = sigmak2;
end

param.Betak = Betak;
param.Sigmak2 = Sigmak2;
%param.Alphak = Alphak;
% else % initialisation aléatoire
%     para.beta =  rand(dimBeta,K);
%     para.sigma = rand(K,1) ;
% end






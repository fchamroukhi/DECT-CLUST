function param = Initialize_SRM(V, B, Y, K, mixingOption, init_kmeans, try_algo)
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

if strcmp(mixingOption, 'softmax')
    %   1. Initialize the mixing weights function's parameters
    if try_algo ==1; param.Alphak = zeros(q1,K-1);else, param.Alphak = rand(q1,K-1);end
    % %   2. Initialization of  and  (from one curve)
    %[param.Betak, param.Sigma2k] = init_regressors(B, Y, K, try_EM);
    
else
    %% Gaussian gating
    
    % Gating Net parameters
    sol_km = myKmeans(V,K);
    klas = sol_km.klas;
    param.Mus = sol_km.muk;
    
    [n, p] = size(V);
    param.R = zeros(p, p, K);
    param.Theta = zeros(p, p, K);
    param.Omega = zeros(1,K);
    for k=1:K
        Vk = V(klas==k, :);
        nk = size(Vk,1);
        param.Omega(k) = nk/n;
        z = (Vk - ones(nk,1)*param.Mus(k,:));
        param.R(:,:,k)  = z'*z/(m*nk);%eye(p);%cov(Xk);
        %z'*z/nk;
        %     if p<n
        %         Theta(:,:,k) = inv(R(:,:,k));
        %     else
        %param.Theta(:,:,k) = eye(p);
        %     end
        
        %param.R(:,:,k)  = eye(p)*sum(sum(z.*z))/(m*nk);
    end
end

%%
dimBeta = size(B, 2);
% Initialize the regression parameters betak and the variances sigma2k

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
    sigmak2 = (1/(nk*m))*sum(sum((Yk-ones(nk,1)*Ykfit').^2));
    %Sigmak2(k)= sigmak2;
    Sigmak2(:,:,k) = sigmak2;%cov(Yk);
end

param.Betak = Betak;
param.Sigmak2 = Sigmak2;
%param.Alphak = Alphak;
% else % initialisation aléatoire
%     para.beta =  rand(dimBeta,K);
%     para.sigma = rand(K,1) ;
% end






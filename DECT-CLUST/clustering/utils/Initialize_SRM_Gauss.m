function param = Initialize_SRM_Gauss(V, Y, K, mixingOption, init_kmeans, try_algo)
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

[n, p]=size(Y);

%basis = regressionOption.basis;
%X = regressionOption.design;

%%   1. Initialize the mixing weights function's parameters
if strcmp(mixingOption, 'softmax')
    if try_algo ==1; param.Alphak = zeros(size(V,2),K-1);else, param.Alphak = rand(size(V,2),K-1);end
else
    %% Gaussian gating
    
    % Gating Net parameters
    sol_km = myKmeans(V,K);
    klas = sol_km.klas;
    param.Mus = sol_km.muk;
    
    [n, d] = size(V);
    param.R = zeros(d, d, K);
    param.Theta = zeros(d, d, K);
    param.Alphak = zeros(1,K);
    for k=1:K
        Vk = V(klas==k, :);
        nk = size(Vk,1);
        param.Alphak(k) = nk/n;
        z = (Vk - ones(nk,1)*param.Mus(k,:));
        %
        lambda = 0.05;
        param.R(:,:,k) = lambda*(z'*z/nk);
        param.Theta(:,:,k) = inv(param.R(:,:,k));
    end
end

%%2. Initialize the regression parameters betak and the variances sigma2k

Muk   = zeros(K, p);
Sigmak2  = zeros(p,p, K);
%Sigmak2  = zeros(m, m, K);

% 2. betak and sigma2k
%init_kmeans=1;
if init_kmeans
    max_iter_kmeans = 400;
    n_tries_kmeans = 3;
    verbose_kmeans = 0;
    sol_kmeans = myKmeans(Y, K,n_tries_kmeans, max_iter_kmeans, verbose_kmeans);
    klas = sol_kmeans.klas;
    Muk = sol_kmeans.muk;
    for k=1:K
        Sigmak2(:,:,k) = cov(Y(klas==k,:)); %diag(diag(cov(Y(klas==k,:))));
    end
else % random partition
    klas = zeros(n,1);
    for i=1:n; klas(i) = randperm(K,1); end
    
    for k=1:K
        Yk = Y(klas==k, :);
        Muk(k,:) = mean(Yk);
        Sigmak2(:,:,k) = cov(Yk);
    end
end


param.Muk = Muk;
param.Sigmak2 = Sigmak2;


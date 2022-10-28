function param = Initialize_SRM(V, B, Y, K, mixingOption, init_kmeans, try_algo, lambda)
%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%%%%%%%%%%%%%%%%%% FC

[n, m]=size(Y);
q1 = size(V,2);


%% 1. Initialize the regression parameters betak and the variances sigma2k
dimBeta = size(B, 2);

Betak   = zeros(dimBeta,K);
Sigmak2  = zeros(K,1);
%Sigmak2  = zeros(m, m, K);

% 2. betak and sigma2k
if init_kmeans
    max_iter_kmeans = 500;  % 300
    n_tries_kmeans = 1;
    verbose_kmeans = 0;
    sol_kmeans = myKmeans(V, K, n_tries_kmeans, max_iter_kmeans, verbose_kmeans);     % spatial !
%     sol_kmeans = myKmeans(Y, K, n_tries_kmeans, max_iter_kmeans, verbose_kmeans);   % curves !
    klas = sol_kmeans.klas;
%     klas = kmeans(Y, K, 'MaxIter', max_iter_kmeans);
else % random partition
    klas = zeros(n,1);
    for i=1:n; klas(i) = randperm(K,1); end
end

for k=1:K
    Yk = Y(klas==k ,:); %if kmeans
    [nk, m] = size(Yk);
    
    Btild = repmat(B,nk,1);
        
    Yktild = reshape(Yk',[],1);
    betak =  (Btild'*Btild)\(Btild'*Yktild);
    Betak(:,k) = betak;    
    Ykfit =  B*betak;
    sigmak2 = sum(sum((Yk-ones(nk,1)*Ykfit').^2))/(m*nk);
    Sigmak2(k) = sigmak2;%cov(Yk);
end

param.Betak = Betak;
param.Sigmak2 = Sigmak2;

%%   2. Initialize the mixing weights function's parameters
if strcmp(mixingOption, 'softmax')
    if try_algo ==1; param.Alphak = zeros(q1,K-1);else, param.Alphak = rand(q1,K-1);end
else
    %% Gaussian gating
    
    % Gating Net parameters
    if init_kmeans
        sol_km = sol_kmeans;
    else
        max_iter_kmeans = 300;
        n_tries_kmeans = 2;
        verbose_kmeans = 0;
        sol_km = myKmeans(V, K , n_tries_kmeans, max_iter_kmeans, verbose_kmeans);   % kmeans for spatial coordinates !
%         sol_km = myKmeans(Y, K , n_tries_kmeans, max_iter_kmeans, verbose_kmeans);   % kmeans for curves !
    end
%     save(['solInitKmeans_cl',num2str(K),'_',num2str(randi(1000)),'.mat'],'sol_km');
%     save(['solInitKmeansCurves_cl',num2str(K),'_',num2str(randi(1000)),'.mat'],'sol_kmeans');
    klas = sol_km.klas;
    param.Mus = sol_km.muk;
    
    [n, p] = size(V);
    param.R = zeros(p, p, K);
    param.Theta = zeros(p, p, K);
    param.Alphak = zeros(1,K);
    for k=1:K
        Vk = V(klas==k, :);
        nk = size(Vk,1);
        param.Alphak(k) = nk/n;
        %param.Mus(k,:) = mean(Vk);
        
        z = (Vk - ones(nk,1)*param.Mus(k,:));
        
        %param.R(:,:,k) = lambda*diag(diag(z'*z/nk));
        param.R(:,:,k)  = lambda*(z'*z)/(nk);
        param.Theta(:,:,k) = inv(param.R(:,:,k));
    end
end


function [klas, params, Posterior, gmm_density, stored_loglik] = learn_RE_SRM_EM(x, Y, K, p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Robust EM algorithm for Random Effects Polynomial Regression Mixture Model
%
%
%
%
%
% by Faicel Chamroukhi, December 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% warning off all

[n, m] = size(Y);
Y_in  = Y;

% Construction of the desing matrix
X = designmatrix_Poly_Reg(x,p);%[m x (p+1)]


%n regularly sampled curves
Xstack = repmat(X,n,1);% desing matrix [(n*m) x (p+1)]


Epsilon = 1e-6;

%------ Step 1 initialization ----- %

% -----------(28) ---------------- %
Y_stack_tmp = repmat(Y',n,[]); % [mxn x n]
Y_stack     = (reshape((Y_stack_tmp(:))',m,[]))'; %[nxn x m];
dij         = sum((Y_stack - repmat(Y,n,[])).^2, 2);
dmin        = min(dij(dij>0));
Q           = dmin;

%%%%
Ytild = reshape(Y',[],1); % []
%%%

%Initialize the mixing proportins
Pik = 1/K*ones(K,1);

% Initialize the regression parameters and the variances
Betak   = zeros(p+1,K);
%d = p+1;
Sigmak2  = zeros(K,1);
for k=1:K
    % ------- step 2  (27)  ------- %
    %betak  = inv(Phi'*Phi + 1e-4*eye(p+1))*Phi'*Y_in(k,:)';
    betak  = (X'*X)\(X'*Y_in(k,:)');
    Betak(:,k) = betak;
    muk = X*betak;
    %Dk = sum((reshape(X,n,m) - reshape(muk',n,m)).^2, 2);
    Dk = sum((Y_in - ones(n,1)*muk').^2, 2);
    Dk = sort(Dk);
    Sigmak2(k)=  Dk(ceil(sqrt(K)));%sum(Y_in(k,:)' - muk);%Dk(ceil(sqrt(K)));%.001;%%;%
    % --------------------------- %
end


% param = initialize_MixReg(Y, K , S, init_kmeans, try_EM);
% init_kmeans = 0;
% try_EM = 1;
% param = initialize_MixReg(Y, K , Xstack, init_kmeans, try_EM);
% 
% Pik = param.Pi_k;
% Betak = param.beta_k;
% Sigmak2 = param.sigma_k;
% RE
Ksi_k = rand(K,1);
Sigmak2 = rand(K,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           %
% main EM-RE-PRM loop       %
%                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stored_loglik  = []; % to store the maximized penalized log-likelihood criterion
loglik_old = -inf;
iter      = 1; % iteration number
MaxIter   = 100;
converged = 0;

%RE
wXbk = zeros(n*m,1);
Lambda_ik_tild= zeros(n*m,1);
while(iter<= MaxIter)% && ~converged)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           %
    %       E-Step              %
    %                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % E-Step of Gaffney
    % [N,D] = size(Y);
    % [P,K,D] = size(Mu);
    %
    % Pik = zeros(n,K,D);
    % Beta_ikd = zeros(P,n,K,D);
    % Vikd = zeros(P,P,n,K,D);
    % for k=1:K
    %   for d=1:D
    %     Rinv = inv(R(:,:,k,d));
    %     for i=1:n
    %       indx = Seq(i):Seq(i+1)-1;
    %       n_i = length(indx);
    %       V_y = X(indx,:)*R(:,:,k,d)*X(indx,:)' + Sigma(i,d)*eye(n_i);
    %       Pik(i,k,d) = mvnormpdf(Y(indx,d)',X(indx,:)*Mu(:,k,d),V_y);
    %
    %       A = 1/Sigma(i,d)*X(indx,:)'*X(indx,:) + Rinv;
    %       c = 1/Sigma(i,d)*X(indx,:)'*Y(indx,d) + Rinv*Mu(:,k,d);
    %       Beta_ikd(:,i,k,d) = A\c;
    %
    %       Vikd(:,:,i,k,d) = inv(1/Sigma(i,d)*X(indx,:)'*X(indx,:) + Rinv);
    %     end
    %   end
    % end
    % Pik = prod(Pik,3);  % scaling problems?
    % Pik = Pik .* (ones(n,1)*Alpha');
    %PikFik = zeros(n, K);
    %log_Pik_fk_Xi = zeros(n,K);
    log_Pik_Fik = zeros(n,K);
    % RE
    D       = p+1;
    Ki_ik   = zeros(n,K);
    lambda_ik = zeros(n,K);
    B_ik = zeros(D,n,K);
    %
    X_i = X;
    m_i=m;
    for k=1:K
        pik = Pik(k);
        betak = Betak(:,k);
        sigmak2 = Sigmak2(k);       
        %if heteroscedasticity==1,sigmak  = param.sigma_k(k)  ;%variance %[1] ou [1 x K]
        %else sigmak=param.sigmak;end
        %
        ksi_k = Ksi_k(k);
        %
        Muk = reshape(Xstack*betak, m, n)'; %[n*m];
        
        Sigmaki = ksi_k*(X_i*X_i') + sigmak2*eye(m_i);%[m x m]
         %fik  
        
%          plot(reshape(Xstack*betak, m, n))
%          
%          pause
%          plot(Y)
%          pause
%          plot(Muk)
%          pause

        z =((Y-Muk)/(Sigmaki)).*(Y-Muk);
        mahalanobis = sum(z,2);
        
        log_fk_Xi = - (m_i/2)*log(2*pi) - 0.5*logdet(Sigmaki) - 0.5*mahalanobis;
        log_Pik_Fik(:,k) = log(pik) + log_fk_Xi;% [n x K]
        %
Ki_ik(:,k) =  ones(n,1)*...
    trace(ksi_k*(eye(p+1) - ksi_k*X_i'/(Sigmaki)*X_i));%repmat sur i
               
lambda_ik(:,k) =  ones(n,1)*trace(ksi_k*X_i*(eye(p+1)-ksi_k*X_i'/Sigmaki*X_i)*X_i');%repmat sur i)
    
    B_ik(:,:,k) = ksi_k*(X_i'/Sigmaki)*(Y-Muk)'; %[Dxn xK]
    end
    
    Pik_Fik = exp(log_Pik_Fik);
    %Posterior = PikFik./(sum(PikFik,2)*ones(1,K));
    log_Prosterior  = log_normalize(log_Pik_Fik);
    Posterior       = exp(log_Prosterior);
    Tauik           = Posterior;
    
    %     % print the value of the optimized criterion
    %     %pen_loglik = sum(log(sum(PikFik,2)),1 ) + Beta*n*sum(Alphak.*log(Alphak));
    %     %pen_loglik = (sum(log_Prosterior(:) .* Posterior(:)) - sum(log_Prosterior(:) .* log_Pik_Fik(:)))+ Beta*n*sum(Alphak.*log(Alphak));
    loglik = sum(logsumexp(log_Pik_Fik,2),1);% + Beta*n*sum(Alphak.*log(Alphak));
    fprintf(1,'EM for RE-PRM Iteration : %d | log-likelihood: %f \n',iter, loglik);
    
    
    %%%%%%%%%%
    % compute the value of the optimized criterion J (12) %
    %pen_loglik = sum(log(sum(PikFik,2)),1 ) + Beta*n*sum(Alphak.*log(Alphak));
    %pen_loglik = sum(logsumexp(log_Pik_Fik,2),1) + Beta*n*sum(Alphak.*log(Alphak));
    %pen_loglik = (sum(log_Prosterior(:) .* Posterior(:)) - sum(log_Prosterior(:) .* log_Pik_Fik(:)))+ Beta*n*sum(Alphak.*log(Alphak));
    stored_loglik = [stored_loglik loglik];
    %     fprintf(1,'EM Iteration : %d  | number of clusters K : %d | penalized loglikelihood: %f \n',iter, K, pen_loglik);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           %
    %       M-Step              %
    %                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     close,
%     plot(Tauik)
%     pause
%     close
%     Sigmak2
    for k=1:K
        tauik = Tauik(:,k);
        % update the mixing proportions
        pik = sum(tauik)/n;%alpha_k^EM
        Pik(k) = pik;
        %
        temp =  repmat(tauik,1,m);% [m x n]
        Wk = reshape(temp',[],1);%cluster_weights(:)% [mn x 1]
        % meme chose
        % temp =  repmat(tauik,1,m)';% [m x n]
        % cluster_weights = cluster_weights(:);
        wYk = sqrt(Wk).*Ytild; % fuzzy cluster k
        wXk = sqrt(Wk*ones(1,p+1)).*Xstack;%[(n*m)*(p+1)]
        %% RE
        bk = B_ik(:,:,k);
        % wXbk = sqrt(Wk*ones(1,p+1)).*(Xstack*bk');%[(n*m)*1]
        for i=1:n
            wXbk((i-1)*m+1:i*m) = (sqrt(Wk((i-1)*m+1:i*m)*ones(1,p+1).*Xstack((i-1)*m+1:i*m,:))*bk(:,i));%[(n*m)*1]
        end
        
        % update the regression coefficients
        betak  = (wXk'*wXk)\(wXk'*(wYk - wXbk)); % + 1e-4*eye(p+1)
        Betak(:,k) = betak;
        
  
        %update the kxi_k's
%         Ksi_k(k)= sum (tauik.*(sum(bk.*bk,1)' + Ki_ik(:,k)))/(D*sum(tauik));
        Ksi_k(k)= sum (tauik.*(sum(bk.*bk,1)' + Ki_ik(:,k)))/(D*sum(Wk));
        
        beta_k = repmat(betak,1,n); % ?
         Y_c = Y'- X_i*(beta_k + bk);
         S_Sigma(k) = tauik'*(sum(Y_c.*Y_c)'+ lambda_ik(:,k));
         
         sigmak2 = sum(tauik'*(sum(Y_c.*Y_c)'+ lambda_ik(:,k)))/sum(Wk);

        %update the variances
        %sigmak2
        %lambda_ik_tild = repmat(lambda_ik(:,k)',m,[]);
        %Lambda_ik_tild(:,k) = lambda_ik_tild(:);
        %sigmak2 = sum((wYk - wXk*betak - wXbk).^2 + Lambda_ik_tild(:,k))/sum(Wk);%(n*m);%*sum(Wk)
        Sigmak2(k) = sigmak2;
        %%
    end
    %Sigmak2=sum(S_Sigma)/(n*m) * ones(K,1);

    % -----------step 11 ---------------- %
    % test of convergence
    
    %     distBetak = sqrt(sum((Betak - BetakOld).^2, 2));
    %     if (max(distBetak) < Epsilon || abs((loglik - loglik_old)/loglik_old)<Epsilon);
    if abs((loglik - loglik_old)/loglik_old)<Epsilon;
        converged = 1;
    end
    loglik_old = loglik;
    
    iter=iter+1;
    
end% en of the Robust EM loop

[~, klas] = max (Posterior,[],2);

gmm_density    = sum(Pik_Fik,2);
params.Pik     = Pik;
params.Betak   = Betak;
params.Muk     = X*Betak;
params.Sigmak2 = Sigmak2;
end




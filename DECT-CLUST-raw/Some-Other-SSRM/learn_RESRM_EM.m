function [mixModel, mixStats] = learn_RESRM_EM(data, K, regressionOption, nbr_EM_runs)
%klas, mixModel, Posterior, model, stored_lglik] = learn_SRM_EM(V, x, Y, K, spline_order, regression_type, nknots)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[mixModel, mixStats] = robust_em_PSRM(x, Y, spline_order, 'spline', nknots);
% EM algorithm for Spatial (Polynomial, (B)-Spline) Regression Mixture Model (SRM) with Random
% Effects
%
%
%
%
%
% FC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off all

V = data.spatialcoord;
% Curves
X = data.abscissas;
Y = data.ordinates;
[n, m] = size(Y);

% Constructing the design matrices for the mixing weights and for the regressor components
[Xw, B] = designSRM(V, X, regressionOption);
%Xw = repmat(Xw,n,1);
% q1 = size(Xw, 2);

dimBeta = size(B, 2);
%n regularly sampled curves
Bstack = repmat(B,n,1);% desing matrix [(n*m) x (dimBeta)]
Ytild = reshape(Y',[],1); % [(n*m) x 1]

%
best_loglik = -inf;
EM_run = 1;
while (EM_run <= nbr_EM_runs)
    if nbr_EM_runs>1,fprintf(1, 'EM run n°  %d \n ',EM_run); end
    EM_run = EM_run + 1;
    %% EM Initializaiton
    
    init_kmeans = 0;
    mixModel = Initialize_SRM(Xw, B, Y, K, init_kmeans, EM_run);
    
    Alphak = mixModel.Alphak;
    Betak = mixModel.Betak;
    Sigmak2 = mixModel.Sigmak2;
    piik = logit_model(Alphak,Xw);
    
    % RE
Ksi_k = rand(K,1);
%Sigmak2 = rand(K,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           %
    % main  EM loop %
    %                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    threshold = 1e-6;% threshold for testing the convergence
    
    stored_lglik  = []; % to store the maximized log-likelihood criterion
    loglik_old = -inf;
    iter      = 0; % iteration number
    MaxIter   = 1500;
    converged = 0;
    
    %RE
wXbk = zeros(n*m,1);
Lambda_ik_tild= zeros(n*m,1);
    while(iter<=MaxIter && ~converged)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           %
        %       E-Step              %
        %                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        log_fk_yij = zeros(n*m,K);
        log_Pik_fk_Yi = zeros(n,K);
        
        
            % RE
    D       = dimBeta;
    Ki_ik   = zeros(n,K);
    lambda_ik = zeros(n,K);
    B_ik = zeros(D,n,K);
    %
    X_i = B;
    m_i=m;
        for k=1:K
            %alphak = Alphak(k);
            betak  = Betak(:,k);
            sigmak2 = Sigmak2(k);
% RE SRM
ksi_k = Ksi_k(k);
%
Muk = reshape(Bstack*betak, m, n)'; %[n*m];
Sigmaki = ksi_k*(X_i*X_i') + sigmak2*eye(m_i);%[m x m]

%             %fik = normpdf(X,muk,sigmak); %Gaussian density
%             z=((Ytild-Bstack*betak).^2)/sigmak2;
%             log_fk_yij(:,k) = - 0.5*(log(2*pi)+log(sigmak2)) - 0.5*z;  %[nxm x 1]
%             % log-lik for the n_k curves of cluster k
%             log_fk_Yi =  sum(reshape(log_fk_yij(:,k),m, n),1); % [n x m]:  sum over j=1,...,m: fk_Xi = prod_j sum_k pi_{jk} N(x_{ij},mu_{k},s_{k))
%             
%             %[logp] = logmvnpdf(X,muk,sigmak)
%             
%             

        %
        z =((Y-Muk)/(Sigmaki)).*(Y-Muk);
        mahalanobis = sum(z,2);
        
        log_fk_Yi = - (m_i/2)*log(2*pi) - 0.5*logdet(Sigmaki) - 0.5*mahalanobis;
        %log_Pik_Fik(:,k) = log(pik) + log_fk_Yi;% [n x K]
        log_Pik_fk_Yi(:,k) = log(piik(:,k)) + log_fk_Yi;% [nxK]

% RE SRM
Ki_ik(:,k) =  ones(n,1)*trace(ksi_k*(eye(dimBeta) - ksi_k*X_i'/(Sigmaki)*X_i));%repmat sur i            
lambda_ik(:,k) =  ones(n,1)*trace(ksi_k*X_i*(eye(dimBeta)-ksi_k*X_i'/Sigmaki*X_i)*X_i');%repmat sur i)
B_ik(:,:,k) = ksi_k*(X_i'/Sigmaki)*(Y-Muk)'; %[Dxn xK]
    
        end
        log_Posterior = log_normalize(log_Pik_fk_Yi);
        Posterior = exp(log_normalize(log_Posterior));
        Tauik = Posterior;
        
        %  print the value of the optimized log-likelihood criterion
        loglik = sum(logsumexp(log_Pik_fk_Yi,2),1);
        fprintf(1,'EM Iteration RE-SRM : %d  log-likelihood: %f \n',iter, loglik);
        %loglik = sum(log(sum(PikFik,2))) + softmax.reg_irls;% + regEM;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           %
        %       M-Step              %
        %                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k=1:K
            tauik = Tauik(:,k);
            % update of the regression coefficients
            tmp =  repmat(tauik,1,m);% [m x n]
            Wk = reshape(tmp',[],1);%cluster_weights(:)% [mn x 1]
            % same thing
            % temp =  repmat(tauik,1,m)';% [m x n]
            % cluster_weights = cluster_weights(:);
            wYk = sqrt(Wk).*Ytild; % fuzzy cluster k
            wBk = sqrt(Wk*ones(1,dimBeta)).*Bstack;%[(n*m)*(M+nknots)]
            % maximization w.r.t betak: Weighted least squares
            %%betak  = inv(phik'*phik + 1e-4*eye(dimBeta))*phik'*Yk;
            %betak  = (wBk'*wBk)\(wBk'*wYk);
            %Betak(:,k) = betak;
            
            % update the variance
            %sigmak2 = sum((wYk - wBk*betak).^2)/sum(Wk);%(m*sum(tauik));
            % % sigmak2 = (1-gama)*sigmak2 + gama*Q;
            %Sigmak2(k) = sigmak2;
       %wYk = sqrt(Wk).*Ytild; % fuzzy cluster k
        wXk = wBk;%sqrt(Wk*ones(1,p+1)).*Xstack;%[(n*m)*(p+1)]    
        %% RE
        bk = B_ik(:,:,k);
        % wXbk = sqrt(Wk*ones(1,p+1)).*(Xstack*bk');%[(n*m)*1]
        for i=1:n
            wXbk((i-1)*m+1:i*m) = (sqrt(Wk((i-1)*m+1:i*m)*ones(1,dimBeta).*Bstack((i-1)*m+1:i*m,:))*bk(:,i));%[(n*m)*1]
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
         
         sigmak2 = ...%sum(tauik'*(sum(Y_c.*Y_c)'+ lambda_ik(:,k)))
         S_Sigma(k)/sum(Wk);
         Sigmak2(k) = sigmak2;
        end
        % update the mixing proportions : Alphak
        
%         plot(B)
%         e
        %%  IRLS : Iteratively Reweighted Least Squares
        softmax = IRLS(Xw, Tauik, Alphak);
        Alphak = softmax.W;
        piik = softmax.piik;
        %% End of EM
        iter=iter+1;
        
        % test of convergence
        if (abs(loglik - loglik_old) <=threshold)
            converged = 1;
        end
        % store the loglik values
        
        stored_lglik = [stored_lglik loglik];
        
        loglik_old = loglik;
    end% en of the EM loop
    
    mixdensity   = sum(exp(log_Pik_fk_Yi),2);
    mixModel.Alphak = Alphak;
    mixModel.Betak  = Betak;
    mixModel.Sigmak2 = Sigmak2;
    %
    mixStats.mixingprobs = piik;
    mixStats.posterior_prob = Tauik;
    mixStats.Muk    = B*Betak;
    mixStats.mixdensity = mixdensity;
    mixStats.stored_loglik = stored_lglik;
    
    if loglik > best_loglik
        best_mixModel = mixModel;
        mixStats.loglik = loglik;
        best_mixStats = mixStats;
        best_loglik = loglik;
    end
    
end % En of EM runs
mixModel = best_mixModel;
mixStats = best_mixStats;
[~, mixStats.klas] = max(mixStats.posterior_prob,[],2);

%
if nbr_EM_runs>1;   fprintf(1,'best loglik:  %f\n',mixStats.loglik); end

end




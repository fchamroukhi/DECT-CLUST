function [mixModel, mixStats] = learn_SRM_EM(data, K, mixingOption, regressionOption, nbr_EM_runs)
%klas, mixModel, Posterior, model, stored_lglik] = learn_SRM_EM(V, x, Y, K, spline_order, regression_type, nknots)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[mixModel, mixStats] = robust_em_PSRM(x, Y, spline_order, 'spline', nknots);
% EM algorithm for Spatial (Polynomial, (B)-Spline) Regression Mixture Model (SRM)
% M: Spline order
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
[Xw, B] = designSRM(V, X, mixingOption, regressionOption);
%Xw = repmat(Xw,n,1);
% q1 = size(Xw, 2);

dimBeta = size(B, 2);
%n regularly sampled curves

Bstack = repmat(B,n,1);% desing matrix [(n*m) x (dimBeta)]

% Bstack = [];
% for i=1:n
%     Bstack = [Bstack; [ones(m, 1)*V(i,:) B]];
% end
% dimBeta = size(Bstack, 2);
% B = Bstack(1:m,:);

Ytild = reshape(Y',[],1); % [(n*m) x 1]

%
best_loglik = -inf;
EM_run = 1;
softmax.reg_irls = 0;

while (EM_run <= nbr_EM_runs)
    if nbr_EM_runs>1,fprintf(1, 'EM run n°  %d \n ',EM_run); end
    EM_run = EM_run + 1;
    %% EM Initializaiton
    
    init_kmeans = 0;
    mixModel = Initialize_SRM(Xw, B, Y, K, mixingOption, init_kmeans, EM_run);
    
    if strcmp(mixingOption,'softmax')
        Alphak = mixModel.Alphak;
        piik = logit_model(Alphak,Xw);
    else%'Gaussian'
        Omega =  mixModel.Omega;
        Mus =  mixModel.Mus;
        R =  mixModel.R;%, R, beta0, Beta, sigma2] = initialize_HDMoGATE(X, y, K);
        Piik = zeros(n, K);
        for k=1:K
            [~, gaussiank] = mvgaussian_pdf(Xw, Mus(k,:), R(:,:,k));
            
            Piik(:,k) = Omega(k)*gaussiank;
        end
        piik = Piik./(sum(Piik,2)*ones(1,K));
    end
    Betak = mixModel.Betak;
    Sigmak2 = mixModel.Sigmak2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           %
    % main  EM loop %
    %                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    threshold = 1e-6;% threshold for testing the convergence
    
    stored_lglik  = []; % to store the maximized log-likelihood criterion
    loglik_old = -inf;
    iter      = 0; % iteration number
    MaxIter   = 300;
    converged = 0;
    while(iter<=MaxIter && ~converged)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           %
        %       E-Step              %
        %                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        log_fk_yij = zeros(n*m,K);
        log_Pik_fk_Yi = zeros(n,K);
        for k=1:K
            %alphak = Alphak(k);
            betak  = Betak(:,k);
                        sigmak2 = Sigmak2(k);
%             sigmak2 = Sigmak2(:,:,k);
            
                        %fik = normpdf(X,muk,sigmak); %Gaussian density
                        z=((Ytild-Bstack*betak).^2)/sigmak2;
                        log_fk_yij(:,k) = - 0.5*(log(2*pi)+log(sigmak2)) - 0.5*z;  %[nxm x 1]
                        % log-lik for the n_k curves of cluster k
                        log_fk_Yi =  sum(reshape(log_fk_yij(:,k),m, n),1); % [n x m]:  sum over j=1,...,m: fk_Xi = prod_j sum_k pi_{jk} N(x_{ij},mu_{k},s_{k))
            
            %Sigmak2(:,:,k) = (sqrt(Wk).*(Ytild - Bstack*betak))'*(sqrt(Wk).*(Ytild - Bstack*betak))/sum(Wk);
            
            %[logp] = logmvnpdf(X,muk,sigmak)
            %log_fk_Yi = logmvnpdf(Y,ones(n,m)*(B*betak),Sigmak2(:,:,k));
            
            log_Pik_fk_Yi(:,k) = log(piik(:,k)) + log_fk_Yi';% [nxK]
        end
        log_Posterior = log_normalize(log_Pik_fk_Yi);
        Posterior = exp(log_normalize(log_Posterior));
        Tauik = Posterior;
        
        %[~, mixStats.klas] = max(mixStats.posterior_prob,[],2);
        
        %mnrnd(mixStats.posterior_prob)
        
        
        %  print the value of the optimized log-likelihood criterion
        loglik = 1/n*sum(logsumexp(log_Pik_fk_Yi,2),1) ;%+ softmax.reg_irls;
        fprintf(1,'EM Iteration : %d  log-likelihood: %f \n',iter, loglik);
        %loglik = sum(log(sum(PikFik,2))) ;% + regEM;
        
        %[~, mixStats.klas] = max(Posterior,[],2);
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
            %betak  = inv(phik'*phik + 1e-4*eye(dimBeta))*phik'*Yk;
            betak  = (wBk'*wBk)\(wBk'*wYk);
            Betak(:,k) = betak;
            %B*betak;
            % update the variance
            sigmak2 = sum((wYk - wBk*betak).^2)/sum(Wk);%(m*sum(tauik));
            %             % sigmak2 = (1-gama)*sigmak2 + gama*Q;
            %             Sigmak2(k) = sigmak2;
            
            z=(Y - ones(n,m)*(B*betak)).*(sqrt(tauik)*ones(1,m));
            % the Gaussian cov matrices
            % Sigmak2(:,:,k) = cov(Y(mixStats.klas==k,:))%diag(diag((z'*z)/sum(tauik)))
            
%             Sigmak2(:,:,k) = sigmak2*eye(m);%(sqrt(tauik).*(Y - ones(n,m)*(B*betak)))'*(sqrt(tauik).*(Y - ones(n,m)*(B*betak)))/(m*sum(tauik));
            Sigmak2(k) = sigmak2;%(sqrt(tauik).*(Y - ones(n,m)*(B*betak)))'*(sqrt(tauik).*(Y - ones(n,m)*(B*betak)))/(m*sum(tauik));
        end
        if strcmp(mixingOption,'softmax')
            % update the mixing proportions : Alphak
            %%  IRLS : Iteratively Reweighted Least Squares
            softmax = IRLS(Xw, Tauik, Alphak);
            Alphak = softmax.W;
            piik = softmax.piik;
        else
            
            for k = 1:K
                Tauk = Tauik(:,k);
                % Gating Network
                % The mixing proportions
                Omega(k) = sum(Tauk)/n;
                % the Gaussian means
                 Mus(k,:) = sum(Xw.*(Tauk*ones(1,size(Xw,2))),1)/sum(Tauk);
                 z=(Xw-ones(n,1)*Mus(k,:)).*(sqrt(Tauk)*ones(1,size(Xw,2)));
%                 % the Gaussian cov matrices

sk = sum(sum(z.*z))/(m*sum(Tauk));
                R(:,:,k) = (z'*z)/(m*sum(Tauk));
                
                [~, gaussiank] = mvgaussian_pdf(Xw, Mus(k,:),R(:,:,k));% sk*eye(2));%R(:,:,k));
 %               [~, gaussiank] = mvgaussian_pdf(Xw, Mus(k,:),sk*eye(2));%R(:,:,k));
                Piik(:,k) = Omega(k)*gaussiank;
            end
            piik = Piik./(sum(Piik,2)*ones(1,K));
        end
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
        if strcmp(mixingOption,'softmax')
            
            mixModel.Alphak = Alphak;
        else%'Gaussian'
            mixModel.Omega =   Omega;
            mixModel.Mus =  Mus;
            mixModel.R = R;%
        end
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




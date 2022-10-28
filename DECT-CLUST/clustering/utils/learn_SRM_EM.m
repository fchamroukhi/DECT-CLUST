function [mixModel, mixStats] = learn_SRM_EM(data, K, mixingOption, regressionOption, nbr_EM_runs, lambda, init_kmeans)
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


Ytild = reshape(Y',[],1); % [(n*m) x 1]

%
best_loglik = -inf;
EM_run = 1;

while (EM_run <= nbr_EM_runs)
    if nbr_EM_runs>1,fprintf(1, 'EM run n�  %d \n ',EM_run); end
    EM_run = EM_run + 1;
    %% EM Initializaiton
    
    mixModel = Initialize_SRM(Xw, B, Y, K, mixingOption, init_kmeans, EM_run, lambda);
    
    if strcmp(mixingOption,'softmax')
        Alphak = mixModel.Alphak;
        piik = logit_model(Alphak,Xw);
    else%'Gaussian'
        Alphak =  mixModel.Alphak;
        Mus =  mixModel.Mus;
        R =  mixModel.R;
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
%     store_alpha = {};   % ADDED BY SEGO
%     iter_time = -ones(MaxIter,1);   % ADDED BY SEGO
    
    while(iter<=MaxIter && ~converged)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           %
        %       E-Step              %
        %                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        log_fk_yij = zeros(n*m,K);
        log_Pik_fk_Yi = zeros(n,K);
        log_alpha_Phi_vy = zeros(n,K);
        
        for k=1:K
            %Gaussian density
            z=((Ytild-Bstack*Betak(:,k)).^2)/Sigmak2(k);
            log_fk_yij(:,k) = - 0.5*(log(2*pi)+log(Sigmak2(k))) - 0.5*z;  %[nxm x 1]
            % conditional log-lik for for component k
            log_fk_Yi =  sum(reshape(log_fk_yij(:,k),m, n),1)'; % [n x m]:  sum over j=1,...,m: fk_Xi = prod_j sum_k pi_{jk} N(x_{ij},mu_{k},s_{k))
            
            %            sigmak2 = Sigmak2(:,:,k);
            %
            %              Muk = reshape(Bstack*Betak(:,k),n,m);
            %              Z =((Y-Muk)*inv(sigmak2)).*(Y-Muk);
            %              mahalanobis = sum(Z,2);
            %              log_fk_Yi = - (m/2)*log(2*pi) - 0.5*logdet(sigmak2) - 0.5*mahalanobis;
            %             %Sigmak2(:,:,k) = (sqrt(Wk).*(Ytild - Bstack*betak))'*(sqrt(Wk).*(Ytild - Bstack*betak))/sum(Wk);
            
            if strcmp(mixingOption, 'gaussian')
                % Gating network conditional density
                [log_Phi_V, ~] = mvgaussian_pdf(Xw, Mus(k,:), R(:,:,k));%, 'diagonal');
                if isnan(log_Phi_V)
                    error("mvgaussian_pdf returns NaN.");
                    % once because R(:,:,k) = [0 0 0;0 0 0;0 0 eps]
                    % because sum(Tauik(:,k)) very small
                    % should remove this cluster from every variables?
                end
                % weighted joint loglik for component k
                log_alpha_Phi_vy(:,k) = log(Alphak(k))*ones(n,1) + log_Phi_V + log_fk_Yi;
                
                %log_Pik_fk_Yi(:,k) = log(piik(:,k)) + log_fk_Yi;% [nxK]
                
            else %'softmax'
                % weighted conditional loglik for component k
                log_Pik_fk_Yi(:,k) = log(piik(:,k)) + log_fk_Yi;% [nxK]
            end
        end
        
        if strcmp(mixingOption, 'softmax')
            
            % log_Posterior = log_normalize(log_Pik_fk_Yi);
            log_sum_Pik_fk_Yi = logsumexp(log_Pik_fk_Yi,2);
            log_Posterior = log_Pik_fk_Yi - log_sum_Pik_fk_Yi*ones(1,K);
            
            loglik = 1/n*sum(log_sum_Pik_fk_Yi,1);
            %loglik = 1/n*sum(logsumexp(log_Pik_fk_Yi,2),1);
            
        else
            % log_Posterior = log_normalize(log_alpha_Phi_vy);
            % loglik = 1/m*sum(logsumexp(log_alpha_Phi_vy,2),1);
            
            %log_alpha_Phi_vy = log_normalize(log_alpha_Phi_vy);
            logsum_alpha_Phi_vy = logsumexp(log_alpha_Phi_vy,2);
            log_Posterior = log_alpha_Phi_vy - logsum_alpha_Phi_vy*ones(1,K);
            
            %log_Posterior = log_normalize(log_Pik_fk_Yi);
            
            loglik = 1/n*sum(logsum_alpha_Phi_vy,1);            
        end
        Posterior = exp(log_Posterior);
        %Posterior = exp(log_normalize(log_Posterior));
        Posterior = Posterior./(sum(Posterior,2)*ones(1,K));
        Tauik = Posterior;
        
        % if SEM
        % if iter>10
        % Tauik = mnrnd(1, Tauik);
        % %[Tauik Posterior(1:20,:) sum(Tauik(1:20,:),2)]
        % end
        
        %  print the value of the optimized log-likelihood criterion
%         fprintf(1,'EM Iteration : %d  log-likelihood: %f \n',iter, loglik);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           %
        %       M-Step              %
        %                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k=1:K
            % update of the regression coefficients
            tmp =  repmat(Tauik(:,k),1,m);% [m x n]
            Wk = reshape(tmp',[],1);%cluster_weights(:)% [mn x 1]
            sqrtWk = sqrt(Wk);
            wYk = sqrtWk.*Ytild; % fuzzy cluster k
            %wBk = (sqrtWk*ones(1,dimBeta)).*Bstack;%[(n*m)*(M+nknots)]
            wBk = repmat(sqrtWk,1,dimBeta).*Bstack;%[(n*m)*(M+nknots)]
            % maximization w.r.t betak: Weighted least squares
            betak  = (wBk'*wBk)\(wBk'*wYk);
            Betak(:,k) = betak;
            % update the variance
            sigmak2 = sum((wYk - wBk*betak).^2)/sum(Wk);%(m*sum(tauik));%
            Sigmak2(k) = sigmak2;
        end
        if strcmp(mixingOption,'softmax')
            % update the mixing proportions : Alphak
            %%  IRLS : Iteratively Reweighted Least Squares
            softmax = IRLS(Xw, lambda, Tauik, Alphak);
            Alphak = softmax.W;
            piik = softmax.piik;
        else
            d=size(Xw,2);
            %sk=0;
            for k = 1:K
                % Gaussian mixture gating function
                % The mixing proportions
                nk = sum(Tauik(:,k));
                Alphak(k) = nk/n;
                % the Gaussian means
%                 Mus(k,:) = sum(Xw.*(Tauik(:,k)*ones(1,d)),1)/nk;
                Mus(k,:) = sum(Xw.*(Tauik(:,k)*ones(1,d)),1);
                if nk ~= 0
                    Mus(k,:) = Mus(k,:)/nk;
                end
                Z=(Xw-ones(n,1)*Mus(k,:)).*(sqrt(Tauik(:,k))*ones(1,d));
                
                % the Gaussian cov matrices
%                 R(:,:,k) =  lambda *  (Z'*Z)/(nk);%lambda*diag(diag(Z'*Z/nk));
                R(:,:,k) =  lambda *  (Z'*Z);
                if nk ~= 0
                    R(:,:,k) = R(:,:,k)/nk;
                end
                
            end
        end        
        %% End of EM
        iter=iter+1;
        
        % test of convergence
        if (abs(loglik - loglik_old) <=threshold); converged = 1;end
        % store the loglik values
        
        stored_lglik = [stored_lglik loglik];
        
        loglik_old = loglik;        
        
%         store_alpha{iter} = Alphak;   % ADDED BY SEGO
%         iter_time(iter) = toc;          % ADDED BY SEGO
    end% en of the EM loop
    
    fprintf(1,'End at EM Iteration : %d  log-likelihood: %f \n',iter, loglik);
    
    if strcmp(mixingOption,'softmax')
        
        mixModel.Alphak = Alphak;
        mixStats.mixingprobs = piik;
        
        mixdensity   = sum(exp(log_Pik_fk_Yi),2);
    else%'Gaussian'
        mixModel.Alphak =   Alphak;
        mixModel.Mus =  Mus;
        mixModel.R = R;%
        
        Piik = zeros(n,K);
        for k=1:K%
            [~, gaussiank] = mvgaussian_pdf(Xw, Mus(k,:),R(:,:,k));
            if isnan(gaussiank)
                error("mvgaussian_pdf returns NaN.");
            end
            Piik(:,k) = Alphak(k)*gaussiank;
        end
        piik = Piik./(sum(Piik,2)*ones(1,K));
        mixStats.mixingprobs = piik;
        mixdensity   = sum(exp(log_alpha_Phi_vy),2);
    end
    mixModel.Betak  = Betak;
    mixModel.Sigmak2 = Sigmak2;
    %
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

% save('store_alpha.mat','store_alpha')
% save('iter_time.mat','iter_time')

mixModel = best_mixModel;
mixStats = best_mixStats;
[~, mixStats.klas] = max(mixStats.posterior_prob,[],2);

%
if nbr_EM_runs>1;   fprintf(1,'best loglik:  %f\n',mixStats.loglik); end

end

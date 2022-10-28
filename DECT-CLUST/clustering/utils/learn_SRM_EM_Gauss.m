function [mixModel, mixStats] = learn_SRM_EM_Gauss(data, K, mixingOption, regressionOption, nbr_EM_runs, lambda)
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

% lambda = 0.06;   %5*1e-1;

V = data.spatialcoord;
% Curves
X = data.abscissas;
Y = data.coefficients; % coefficients
[n, m] = size(Y);

% Constructing the design matrices for the mixing weights and for the regressor components
[Xw, ~] =  designSRM(V, X, mixingOption, regressionOption);

best_loglik = -inf;
EM_run = 1;
while (EM_run <= nbr_EM_runs)
    if nbr_EM_runs>1,fprintf(1, 'EM run n°  %d \n ',EM_run); end
    EM_run = EM_run + 1;
    %% EM Initializaiton
    
    init_kmeans = 0;
    mixModel = Initialize_SRM_Gauss(Xw,Y, K, mixingOption, init_kmeans, EM_run);
    
    if strcmp(mixingOption,'softmax')
        Alphak = mixModel.Alphak;
        piik = logit_model(Alphak,Xw);
    else%'Gaussian'
        Alphak =  mixModel.Alphak;
        Mus =  mixModel.Mus;
        R =  mixModel.R;
    end
    Muk = mixModel.Muk;
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
    MaxIter   = 500;
    converged = 0;
    while(iter<=MaxIter && ~converged)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           %
        %       E-Step              %
        %                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        log_Pik_fk_Yi = zeros(n,K);
        log_alpha_Phi_vy = zeros(n,K);
        
        for k=1:K
            
            [log_fk_Yi, ~] = mvgaussian_pdf(Y, Muk(k,:), Sigmak2(:,:,k), 'diagonal');
            if isnan(log_fk_Yi)
                error("mvgaussian_pdf returns NaN.");
            end
            if strcmp(mixingOption, 'gaussian')
                % Gating network conditional density
                [log_Phi_V, ~] = mvgaussian_pdf(Xw, Mus(k,:), R(:,:,k), 'diagonal');
                if isnan(log_Phi_V)
                    error("mvgaussian_pdf returns NaN.");
                end
                % weighted joint loglik for component k
                log_alpha_Phi_vy(:,k) = log(Alphak(k))*ones(n,1) + log_Phi_V + log_fk_Yi;
            else %'softmax'
                log_Pik_fk_Yi(:,k) = log(piik(:,k)) + log_fk_Yi;% [nxK]
            end
        end
        
        if strcmp(mixingOption, 'softmax')
            % log_Posterior = log_normalize(log_Pik_fk_Yi);
            log_sum_Pik_fk_Yi = logsumexp(log_Pik_fk_Yi,2);
            log_Posterior = log_Pik_fk_Yi - log_sum_Pik_fk_Yi*ones(1,K);
            
            loglik = 1/n*sum(log_sum_Pik_fk_Yi,1);
            %loglik = 1/m*sum(logsumexp(log_Pik_fk_Yi,2),1);
        else % normalized Gaussian
            % log_Posterior = log_normalize(log_alpha_Phi_vy);
            % loglik = 1/m*sum(logsumexp(log_alpha_Phi_vy,2),1);
            
            %log_alpha_Phi_vy = log_normalize(log_alpha_Phi_vy);
            logsum_alpha_Phi_vy = logsumexp(log_alpha_Phi_vy,2);
            log_Posterior = log_alpha_Phi_vy - logsum_alpha_Phi_vy*ones(1,K);
            loglik = 1/n*sum(logsum_alpha_Phi_vy,1);
        end
        
        Posterior = exp(log_Posterior);
        Tauik = Posterior./(sum(Posterior,2)*ones(1,K)); % normalize

%         fprintf(1,'EM Iteration : %d  log-likelihood: %f \n',iter, loglik);
        %loglik = sum(log(sum(PikFik,2))) ;% + regEM;
        
        %[~, mixStats.klas] = max(Posterior,[],2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           %
        %       M-Step              %
        %                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k=1:K
            tauik = Tauik(:,k);
            % Gating Network
            % The mixing proportions
            nk = sum(tauik);
            % the Gaussian means
            Muk(k,:) = sum(Y.*(tauik*ones(1,m)),1)/sum(tauik);
            Z = (Y-ones(n,1)*Muk(k,:)).*(sqrt(tauik)*ones(1,m));
            %                 % the Gaussian cov matrices
            Sigmak2(:,:,k) = Z'*Z/nk ;%+ 1e-2*inv(B'*B);%
            
            %sk = sum(sum(Z.^2, 1))/(nk);
            %Sigmak2(:,:,k) = sk*eye(m);%
        end
        if strcmp(mixingOption,'softmax')
            % update the mixing proportions : Alphak
            %%  IRLS : Iteratively Reweighted Least Squares
            softmax = IRLS(Xw, Tauik, Alphak);
            Alphak = softmax.W;
            piik = softmax.piik;
        else
            d=size(Xw,2);
            for k = 1:K
                tauik = Tauik(:,k);
                % Gating Network
                % The mixing proportions
                nk = sum(tauik);
                Alphak(k) = nk/n;
                % the Gaussian means
                Mus(k,:) = sum(Xw.*(tauik*ones(1,d)),1)/sum(tauik);
                Z=(Xw-ones(n,1)*Mus(k,:)).*(sqrt(tauik)*ones(1,d));
                % the Gaussian cov matrices
                R(:,:,k) =  lambda *  (Z'*Z/nk);%lambda*diag(diag(Z'*Z/nk));
                %R(:,:,k) =  ((Z'*Z/(nk)) + lambda *  (Xw'*Xw)/(d*n));
            end
        end
        %% End of EM
        iter=iter+1;
        
        % test of convergence
        if (abs(loglik - loglik_old) <=threshold); converged = 1; end
        % store the loglik values
        
        stored_lglik = [stored_lglik loglik];
        
        loglik_old = loglik;
    end% en of the EM loop
    
    fprintf(1,'End at EM Iteration : %d  log-likelihood: %f \n',iter, loglik);
    
    mixdensity   = sum(exp(log_Pik_fk_Yi),2);
    if strcmp(mixingOption,'softmax')
        
        mixModel.Alphak = Alphak;
        mixStats.mixingprobs = piik;
    else%'Gaussian'
        mixModel.Alphak =   Alphak;
        mixModel.Mus =  Mus;
        mixModel.R = R;%
        
        Piik = zeros(n,K);
        for k=1:K%
            [~, gaussiank] = mvgaussian_pdf(Xw, Mus(k,:),R(:,:,k),'diagonal');
            if isnan(gaussiank)
                error("mvgaussian_pdf returns NaN.");
            end
            Piik(:,k) = Alphak(k)*gaussiank;
        end
        piik = Piik./(sum(Piik,2)*ones(1,K));
        mixStats.mixingprobs = piik;
    end
    mixModel.Muk  = Muk;
    mixModel.Sigmak2 = Sigmak2;
    %
    mixStats.posterior_prob = Tauik;
    %mixStats.Muk    = B*Betak;
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
function [mixModel, mixStats] = learn_SRMF_EM(data, K, mixingOption, regressionOption, nbr_EM_runs)
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
[Xw, B] = designSRMF(V, X, mixingOption, regressionOption);
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

while (EM_run <= nbr_EM_runs)
    if nbr_EM_runs>1,fprintf(1, 'EM run n°  %d \n ',EM_run); end
    EM_run = EM_run + 1;
    %% EM Initializaiton
    
    init_kmeans = 0;
    mixModel = Initialize_SRMF(V, B, Y, K, mixingOption, init_kmeans, EM_run);
    
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
    MaxIter   = 100;
    converged = 0;
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
            %alphak = Alphak(k);
            betak  = Betak(:,k);
            sigmak2 = Sigmak2(k);
            %             sigmak2 = Sigmak2(:,:,k);
            
            %fik = normpdf(X,muk,sigmak); %Gaussian density
            z=((Ytild-Bstack*betak).^2)/sigmak2;
            log_fk_yij(:,k) = - 0.5*(log(2*pi)+log(sigmak2)) - 0.5*z;  %[nxm x 1]
            % log-lik for the n_k curves of cluster k
            log_fk_Yi =  sum(reshape(log_fk_yij(:,k),m, n),1)'; % [n x m]:  sum over j=1,...,m: fk_Xi = prod_j sum_k pi_{jk} N(x_{ij},mu_{k},s_{k))
            
            
            % Muk = reshape(Bstack*betak,n,m);
            %
            % Z =((Y-Muk)*inv(sigmak2)).*(Y-Muk);
            % mahalanobis = sum(Z,2);
            %
            % log_fk_Yi = - (m/2)*log(2*pi) - 0.5*logdet(sigmak2) - 0.5*mahalanobis;
            
            %Sigmak2(:,:,k) = (sqrt(Wk).*(Ytild - Bstack*betak))'*(sqrt(Wk).*(Ytild - Bstack*betak))/sum(Wk);
            
            if strcmp(mixingOption, 'gaussian')
                % Gating network conditional density
                Phi_V = zeros(n,1);
                for i=1:n
                    V_i = Xw{i};
                    %for j=1:lngth(V_i);
                        %[log_Phi_Vij, ~] = mvgaussian_pdf(V_i(j,:), Mus(k,:), R(:,:,k));%, 'diagonal');
                    %end
                    [~, Phi_vij] = mvgaussian_pdf(V_i, Mus(k,:), R(:,:,k));%, 'diagonal');
                    Phi_Vi = sum(Phi_vij,1);
                    %Phi_Vi = sum(mvnpdf(V_i, Mus(k,:), R(:,:,k), 'diagonal'), 1);%, 'diagonal');
                    Phi_V(i) = Phi_Vi;
                end
                log_Phi_V = log(Phi_V);
                % Expert Network conditional density
                %log_Phi_y =  -0.5*log(2*pi) - 0.5*log(sigma2(k)) -0.5*((y - beta0(k)*ones(n,1) - X*Beta(:,k)).^2)/sigma2(k);
                %log_Phi_y =  gaussian_pdf(y, beta0(k)*ones(n,1) + X*Beta(:,k), sigma2(k));
                
                % weighted MoGATE conditional density
                %log_Pik_fk_Yi(:,k) = log(Alphak(k))*ones(n,1) + log_Phi_V + log_fk_Yi;
                log_alpha_Phi_vy(:,k) = log(Alphak(k))*ones(n,1) + log_Phi_V + log_fk_Yi;
                
                %[logp] = logmvnpdf(X,muk,sigmak)
                %log_fk_Yi = logmvnpdf(Y,ones(n,m)*(B*betak),Sigmak2(:,:,k));
            else %'softmax'
                log_Pik_fk_Yi(:,k) = log(piik(:,k)) + log_fk_Yi;% [nxK]
            end
        end
        
        if strcmp(mixingOption, 'softmax')
            
            log_Posterior = log_normalize(log_Pik_fk_Yi);
            Posterior = exp(log_normalize(log_Posterior));
            Tauik = Posterior;
            %
            loglik = 1/n*sum(logsumexp(log_Pik_fk_Yi,2),1);
        else
            log_Posterior = log_normalize(log_alpha_Phi_vy);
            Posterior = exp(log_normalize(log_Posterior));
            Tauik = Posterior;
            %log_alpha_Phi_vy = log_normalize(log_alpha_Phi_vy);
            
            %                 log_sum_alpha_Phivy = logsumexp(log_alpha_Phi_vy,2);
            %         log_Tau = log_alpha_Phi_vy - log_sum_alpha_Phivy*ones(1,K);
            %         Tauik = exp(log_Tau);
            %         Tauik = Tauik./(sum(Tauik,2)*ones(1,K));
            
            %
            %       loglik = 1/m*sum(log_sum_alpha_Phivy);
            loglik = 1/n*sum(logsumexp(log_alpha_Phi_vy,2),1);
        end
        
        %[~, mixStats.klas] = max(mixStats.posterior_prob,[],2);
    
        %  print the value of the optimized log-likelihood criterion
        %loglik = 1/m*loglik + softmax.reg_irls;
        
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
            sigmak2 = sum((wYk - wBk*betak).^2)/sum(Wk);%(m*sum(tauik));%
            %             % sigmak2 = (1-gama)*sigmak2 + gama*Q;
            Sigmak2(k) = sigmak2;
            
            %z=(Y - ones(n,m)*(B*betak)).*(sqrt(tauik)*ones(1,m));
            % the Gaussian cov matrices
            % Sigmak2(:,:,k) = cov(Y(mixStats.klas==k,:))%diag(diag((z'*z)/sum(tauik)))
            
            %             Sigmak2(:,:,k) = sigmak2*eye(m);%(sqrt(tauik).*(Y - ones(n,m)*(B*betak)))'*(sqrt(tauik).*(Y - ones(n,m)*(B*betak)))/(m*sum(tauik));
            %Sigmak2(k) = sigmak2;%
            %              Sigmak2(:,:,k) = (sqrt(tauik).*(Y - ones(n,m)*(B*betak)))'*(sqrt(tauik).*(Y - ones(n,m)*(B*betak)))/(sum(tauik));
            %Sigmak2(:,:,k) = (sqrt(tauik).*(Y - reshape(Bstack*betak,n,m)))'*(sqrt(tauik).*(Y - reshape(Bstack*betak,n,m)))/sum(tauik)
            %Sigmak2(:,:,k) = diag(diag(Sigmak2(:,:,k)));
        end
        if strcmp(mixingOption,'softmax')
            % update the mixing proportions : Alphak
            %%  IRLS : Iteratively Reweighted Least Squares
            softmax = IRLS(Xw, Tauik, Alphak);
            Alphak = softmax.W;
            piik = softmax.piik;
        else
            d=size(Xw,2);
            sk=0;
            for k = 1:K
                tauik = Tauik(:,k);
                % Gating Network
                % The mixing proportions
                nk = sum(tauik);
                Alphak(k) = nk/n;
                % the Gaussian means
%                 Mus(k,:) = sum(Xw.*(tauik*ones(1,d)),1)/sum(tauik);
%                 Z=(Xw-ones(n,1)*Mus(k,:)).*(sqrt(tauik)*ones(1,d));
%                 %                 % the Gaussian cov matrices
%                 %sk = sk+sum(Z.^2, 1)/nk;
%                 sk = sum(sum(Z.^2, 1))/(m*nk);
%                 R(:,:,k) = sk*eye(d) + 1e-6*eye(d);%
%                 %R(:,:,k) = diag(diag((Z'*Z)/(m*nk)));
                
                muk = 0;
                for i=1:n
                    V_i = Xw{i};
                    %muk = muk + Tauik(i,k)*V_i(j,:);  
                    muk = muk + sum(Tauik(i,k)*V_i, 1)/(length(V_i)*nk);
                end
                Mus(k,:) = muk;
                
                sigk= 0;
                
                lambda = 0.1;%0.07;
                for i=1:n
                    V_i = Xw{i};
                    v = length(V_i);
                Z= sqrt(Tauik(i,k))* (V_i- ones(v,1)*Mus(k,:));
                sigk = sigk + lambda * Z'*Z/(length(V_i)*nk);
                %sigk = sigk + sum(sum(Z.^2, 1))/(length(V_i)*nk);
                end
                %R(:,:,k) = sigk*eye(size(V_i, 2));%sigk;
                R(:,:,k) =  sigk;
                
            end

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
        mixStats.mixingprobs = piik;
    else%'Gaussian'
        mixModel.Alphak =   Alphak;
        mixModel.Mus =  Mus;
        mixModel.R = R;%
        
%        Piik = zeros(n,K);
%         for k=1:K%
%             [~, gaussiank] = mvgaussian_pdf(Xw, Mus(k,:),R(:,:,k));
%             Piik(:,k) = Alphak(k)*gaussiank;
%         end
%         piik = Piik./(sum(Piik,2)*ones(1,K));
%         mixStats.mixingprobs = piik;
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
mixModel = best_mixModel;
mixStats = best_mixStats;
[~, mixStats.klas] = max(mixStats.posterior_prob,[],2);

%
if nbr_EM_runs>1;   fprintf(1,'best loglik:  %f\n',mixStats.loglik); end

end



function [log_fxi, fxi] = mvgaussian_pdf(X, mu, sigma, covtype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Clacul d'une densite gaussienne  pour un echantillon donné
% Entrees:
%
%           X    : Tableau d'individus [nxd]
%           mu [1xd]; centre
%           sigma % mat de variance [dxd]
% Sorties:
%
%            fxi : (2*pi*|sigma|)^(-.5*d)*exp{-0.5*(xi - mu)
%            sigma^-1(xi-mu)'}  pour i=1,...n. vecteur de dim n
%            log_fxi:  log(fxi) pour i=1,...n. vecteur de dim n pour
%
%Faicel Chamroukhi
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin>3
    switch covtype
        case 'diagonal'
            %             sigma
            diagelements = diag(sigma);
            %             pause
            InvSigma = diag(1./diagelements);
            %detSigma = prod(diagelements);
            logdetSigma = sum(log(diagelements));
        case 'full'
            %detSigma = det(sigma);
            logdetSigma = logdet(sigma);
            InvSigma = inv(sigma);
        case 'precision'%sigma is the precision matrix
            InvSigma = sigma;
            %detSigma = 1/det(sigma);
            logdetSigma = -logdet(sigma);
        otherwise %full
            %detSigma = det(sigma);
            logdetSigma = logdet(sigma);
            InvSigma = inv(sigma);
    end
else
    try
        %detSigma = det(sigma);
        logdetSigma = logdet(sigma);
        InvSigma = inv(sigma);
    catch
        diagelements = diag(sigma);
        %detSigma = prod(diagelements);
        InvSigma = diag(1./diagelements);
        logdetSigma = sum(log(diagelements));
    end
end
[n, d]=size(X);


z =((X-ones(n, 1)*mu)*InvSigma).*(X-ones(n, 1)*mu);
mahalanobis = sum(z,2);
log_fxi = - (d/2)*log(2*pi) - 0.5*logdetSigma - 0.5*mahalanobis;

% denom = (2*pi)^(d/2)*(detSigma)^(1/2);
% fxi =  exp(-0.5*mahalanobis) /denom;

fxi = exp(log_fxi);
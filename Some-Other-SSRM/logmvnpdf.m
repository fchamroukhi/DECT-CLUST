function [logp] = logmvnpdf(X,muk,sigmak)
% outputs log likelihood array for observations x  where x_n ~ N(mu,Sigma)
% x is NxD, mu is 1xD, Sigma is DxD

[n,d] = size(X);
% X_centered = X-ones(n,1)*muk;
X_centered = X-muk;
% logp = ((2*pi)*(-d/2))+(logdet(sigmak)*(-1/2)) +(-1/2)*sum((X_centered*(sigmak^(-1))).*X_centered,2);

invSigma= (sigmak + 1e-4*eye(size(sigmak)))\eye(size(sigmak));
logp = (-.5*d*log(2*pi)) - (.5 * logdet(sigmak)) - (.5 * sum((X_centered*invSigma).*X_centered,2));

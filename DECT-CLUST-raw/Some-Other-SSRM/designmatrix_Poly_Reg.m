function X = designmatrix_Poly_Reg(x,p)
%
%
%
%
%
%
%
%
%
%
%
%
%
%
% Faicel Chamroukhi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(x,1) == 1;
    x=x'; % en vecteur
end
X=[];
for i = 0:p
    X =[X x.^i];% [1 x x.^2 x.^3 tx.^p;......;...]
end
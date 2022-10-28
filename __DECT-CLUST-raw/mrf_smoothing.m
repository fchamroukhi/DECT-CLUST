function [klas, K] = mrf_smoothing(coord, klas, neighb)

n = length(klas);
K = max(klas);
C = zeros(n, 3);
for i=1:n
    si1 = coord(i,1); si2 = coord(i,2); % in 2d
    %     delta_ij = zeros(n,1);
    %         for j=1:n
    %             sj1 = coord(j,1); sj2 = coord(j,2);
    %             delta_ij(j) = max(abs(si1 - sj1),abs(si2 - sj2));
    %         end
    Sj1 = coord(:,1); Sj2 = coord(:,2);
    delta_ij = max(abs(si1 - Sj1),abs(si2 - Sj2));
    
    Si = find(delta_ij  <= neighb);
    Si (Si == i) = [];
    
    Cs = klas(Si);
    muCs_i = sum(Cs == klas(i))/length(Si); %cardSi = length(Si);
    
    %muCs_i = mean(coord(Cs,:));
    C(i, :) = [1, klas(i) muCs_i];
    %C(i, :) = [1, muCs_i];
end

Zik = repmat(klas,1,K)==(ones(n,1)*[1:K]);
softmax = IRLS(C, Zik);
%W = softmax.W;
piik = softmax.piik;

% % equivalent to
% W = mnrfit(C(:,2:end), klas);
% piik = mnrval(W, C(:,2:end));

[~, smoothed_data]  = max(piik,[],2);
%data = [];
Ksmooth = unique(smoothed_data);
for k=1:length(Ksmooth)
    klas(smoothed_data==Ksmooth(k)) = k;
end
K = length(Ksmooth);
end

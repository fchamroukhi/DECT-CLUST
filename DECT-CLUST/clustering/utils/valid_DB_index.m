function res = valid_DB_index(data,labels,tum_klas_idx)

% adaptation of Davis Bouldin clustering index for tumor fitting 
% Segolene Brivet, segolene.brivet@mail.mcgill.ca, Dec 2021


K = max(labels);
DB = zeros(K,1);

% For tumor merged clusters
tum_ind = labels == tum_klas_idx;
tum_center = mean(data(tum_ind,:),1);  % cluster center
Stum = mean( sqrt(sum( (data(tum_ind,:)-repmat(tum_center,sum(tum_ind),1)).^2, 2)) );  % intra-variability

% For other clusters
for k=1:K
    if k==tum_klas_idx, DB(k) = NaN; continue; end
    
    clust_ind = labels == k;
    if sum(clust_ind) == 0
        DB(k) = NaN;
    else
        clust_center = mean(data(clust_ind,:),1);  % cluster center
        Sk = mean( sqrt(sum( (data(clust_ind,:)-repmat(clust_center,sum(clust_ind),1)).^2, 2)) );  % intra-variability


        DB(k) = (Stum + Sk) / norm(tum_center-clust_center,2);  % clustering score between 2 clusters
    end
end


res = DB;
% res = max(DB);  % Final clustering score




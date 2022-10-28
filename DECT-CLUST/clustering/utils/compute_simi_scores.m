
function [dice_array, jacc_array] = compute_simi_scores(klas, klas_tum, K)

dice_array = zeros(1,K);  % Dice score
jacc_array = zeros(1,K);  % Jaccard = IoU

for cl_id=1:K
    ind = klas==cl_id;
    dice_array(cl_id) = dice(ind,klas_tum);
    jacc_array(cl_id) = jaccard(ind,klas_tum);  % = sum(ind & klas_tum) / sum(ind | klas_tum);
end

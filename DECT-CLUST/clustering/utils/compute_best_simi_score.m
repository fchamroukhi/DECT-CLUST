
function [dice_score, jacc_score] = compute_best_simi_score(klas, klas_tum, match_klas)

ind = false(size(klas));   % logical zeros

for cl_id=match_klas
    ind = ind | klas==cl_id;
end

dice_score = dice(ind,klas_tum);
jacc_score = jaccard(ind,klas_tum);  % = sum(ind & klas_tum) / sum(ind | klas_tum);

end
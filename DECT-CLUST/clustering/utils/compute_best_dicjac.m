
function [best_dice, best_jacc, sel_cl] = compute_best_dicjac(klas_tum, klas, K)


% Note clusters atha covers more than 5% of tumor
int_score = [];
match_klas = [];

clust_ids = unique(klas);
clust_ids = clust_ids(clust_ids>0);

for cl_id = clust_ids'
    curr_klas = (klas == cl_id);
    
    intrs = sum(klas_tum & curr_klas, 'all');
    
    if intrs > sum(klas_tum,'all')*0.05   % if cluster intersects at least 5% of Ground Truth
        int_score = [int_score, intrs];
        match_klas = [match_klas, cl_id];
    end
    
end


% Sort them by biggest coverage
[~, ind] = sort(int_score, 'descend');
match_klas = match_klas(ind);


% Compute Dice score by adding cluster by cluster, and select the best combination
un_cl = false(size(klas));  % logical zeros
best_dice = 0;
sel_cl = [];

for cl_id=match_klas
    sel = un_cl | klas==cl_id;
    
    curr_dice = dice(sel, klas_tum);
    if curr_dice > best_dice
        best_dice = curr_dice;
        sel_cl = [sel_cl, cl_id];
        un_cl = sel;
    end
end
% 
% un_cl = false(size(klas));  % logical zeros
% best_dice = 0;
% fin_cl = [];
% for cl_id=flip(sel_cl)
%     sel = un_cl | klas==cl_id;
%     
%     curr_dice = dice(sel, klas_tum);
%     if curr_dice > best_dice
%         best_dice = curr_dice;
%         fin_cl = [fin_cl, cl_id];
%         un_cl = sel;
%     end
% end
% sel_cl = fin_cl;

best_jacc = jaccard(un_cl, klas_tum);

end

function match_klas = match_tumor_clusters(klas_tum,klas,K)

match_klas = [];
for cl_id=1:K
    curr_klas = (klas == cl_id);
    if sum(klas_tum & curr_klas) > sum(klas_tum)*0.2   % if cluster intersects at least 40% of Ground Truth (before = 10%)
        match_klas = [match_klas, cl_id];
    end
end


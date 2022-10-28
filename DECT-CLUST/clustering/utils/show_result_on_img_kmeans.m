
function show_result_on_img_kmeans(T, img_4D, focus, vars)

nk = length(vars.match_klas);
cmap = hot(8+nk);
cmap = cmap(2:2+nk-1,:);

for i=1:min(vars.max_slic_subplot,length(vars.slic_show))
    
    slic_r = img_4D(:,:,vars.slic_show(i),1);
    slic_g = slic_r; slic_b = slic_r;
    
    subplot(vars.len_sp,4,2*i);
    for cl_id=1:vars.K
        cc = find(T(:,:,i) == cl_id);
        tum_idx = find(cl_id==vars.match_klas);
        if ~isempty(tum_idx)    % shade of red colors for tumor clusters
            slic_r(cc) = cmap(tum_idx,1);
            slic_g(cc) = cmap(tum_idx,2);
            slic_b(cc) = cmap(tum_idx,3);
        else                    % random colors for other clusters
            slic_r(cc) = vars.clr(cl_id,1); slic_g(cc) = vars.clr(cl_id,2); slic_b(cc) = vars.clr(cl_id,3);
        end
    end
    slic_rgb = cat(3,slic_r, slic_g, slic_b);
    imshow(slic_rgb.*focus(:,:,i))
    hold on
    plot_tumor_contour(vars.tumor_contour_list{i}, [vars.rmin, vars.cmin], [0.99,0.99,0.99]);
end

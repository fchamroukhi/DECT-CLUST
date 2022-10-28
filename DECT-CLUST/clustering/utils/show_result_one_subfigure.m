
function show_result_one_subfigure(reco_lbl, vars, i)

nk = length(vars.match_klas);
cmap = hot(8+nk);
cmap = cmap(2:2+nk-1,:);

slic_r = map_to_0_1_wdw(vars.subj_slic(:,:,vars.slic_show(i)),vars.mnx);
slic_g = slic_r; slic_b = slic_r;

for cl_id=1:vars.K
    cc = find(reco_lbl(:,:,vars.slic_show(i)) == cl_id);
    tum_idx = find(cl_id==vars.match_klas);
    if ~isempty(tum_idx)  % shade of red colors
        slic_r(cc) = cmap(tum_idx,1);
        slic_g(cc) = cmap(tum_idx,2);
        slic_b(cc) = cmap(tum_idx,3);
    else
        slic_r(cc) = vars.clr(cl_id,1); slic_g(cc) = vars.clr(cl_id,2); slic_b(cc) = vars.clr(cl_id,3);
    end
end
slic_rgb = cat(3,slic_r, slic_g, slic_b);
imshow(slic_rgb(vars.rmin:vars.rmax,vars.cmin:vars.cmax,:),[]);
hold on
plot_tumor_contour(vars.tumor_contour_list{i}, [vars.rmin, vars.cmin], [0.99,0.99,0.99]);

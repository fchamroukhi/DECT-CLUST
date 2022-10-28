

function fig_slic = plot_classif_results(reco_lbl, tumor_contour_list, subj_slic, slic_show, max_slic_subplot, slic_min_idx, lvl, wdw, rmin, rmax, cmin, cmax)

% subj_slic(subj_slic<(lvl-(wdw/2))) = lvl-(wdw/2);
% subj_slic(subj_slic>(lvl+(wdw/2))) = lvl+(wdw/2);
% slics = subj_slic(rmin:rmax,cmin:cmax,:);
len_sp = min(ceil(length(slic_show)/2),ceil(max_slic_subplot/2));
    
fig_slic = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(slic_show)
    if i > max_slic_subplot, break; end
    
    
    % Tissue enhanced 40kev image
    slic_r = map_to_0_1_wdw(subj_slic(:,:,slic_show(i)),[lvl-wdw/2, lvl+wdw/2]);
    slic_g = slic_r; slic_b = slic_r;
    
    
    % Plot original image
    subplot(len_sp,4,2*i-1);
    imshow(slic_r(rmin:rmax,cmin:cmax),[]), hold on, title("slice " + num2str(slic_show(i)-1+slic_min_idx)) % index showed on 3DSlicer
    plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);

    
    % Plot corresponding results image
    subplot(len_sp,4,2*i);
    
    cc = find(reco_lbl(:,:,slic_show(i)) == 1);         % tumor
    slic_r(cc) = 0.83; slic_g(cc) = 0.14; slic_b(cc) = 0.14;     % red
    
    cc = find(reco_lbl(:,:,slic_show(i)) == 0);         % non-tumor
    slic_r(cc) = 0.8; slic_g(cc) = 0.8; slic_b(cc) = 1; % white-blue
    
%     cc = find(reco_lbl(:,:,slic_show(i)) == -1);        % not classified
%     slic_r(cc) = 0.3; slic_g(cc) = 0.1; slic_b(cc) = 0.4; % brown
    
    slic_rgb = cat(3,slic_r, slic_g, slic_b);
    
    imshow(slic_rgb(rmin:rmax,cmin:cmax,:),[])
    hold on
    plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
end

end

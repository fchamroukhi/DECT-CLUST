
    
%% Ground Truth tumor curve details

for patient_name = ["subject8_tumor","HNSCC2","HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10"]
% for patient_name = ["subject8_tumor","HNSCC2","HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10","HNSCC11","HNSCC12","HNSCC13","HNSCC15","HNSCC17","HNSCC18","HNSCC26"]
    close all; 
    
    load(['../../data/',char(patient_name),'_GT.mat']);
    load(['../../data/',char(patient_name),'.mat']);
    [~,~,slic_min] = ind2sub(size(segm_vol_full),find(segm_vol_full,1,'first')); % lower slice containing a tumor
    [~,~,slic_max] = ind2sub(size(segm_vol_full),find(segm_vol_full,1,'last')); % lower slice containing a tumor
    slic_inds = round((slic_max-slic_min+1)/2)-3 : round((slic_max-slic_min+1)/2)+2; % 6 slices
    %fn_save_pdf = '';  % if don't want to save fig as pdf
    fn_save_pdf = fullfile('GT_curves', char(patient_name));  % to save fig as pdf


    
    %% TUMOR

    % Select tumor in the 3D volume: line 'i' stores the corresponding coordinates and curves
    gr_truth = segm_vol_full(:,:,slic_min + slic_inds -1);
    lin_tum = find(gr_truth);
    [row_obj, col_obj, z_obj] = ind2sub(size(gr_truth),lin_tum);
    
    decay_curves_tum = zeros(length(row_obj),21);
    for r=1:length(row_obj)
        for kev=1:21
            decay_curves_tum(r,kev) = subject{kev}(col_obj(r),row_obj(r),z_obj(r));
        end
    end
    
%     Ytum = (decay_curves_tum - moy*ones(length(decay_curves_tum),1))./stdev;
%     figure;
%     plot(T,decay_curves_tum','color',[1 0 0],'linewidth',0.001);
    

    %% Plot tumor curves
    
    fig_tum = figure; hold on
    plot(40:5:140,decay_curves_tum','color',[1 0 0],'linewidth',0.001);
    xlabel('Energy level','FontWeight','bold');
    ylabel('HU','FontWeight','bold');
    box on;
    title('Curve details for Ground Truth tumor')

    
    if ~isempty(fn_save_pdf)
        set(fig_tum,'Units','Inches');
        pos = get(fig_tum,'Position');
        set(fig_tum,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(fig_tum,[fn_save_pdf,'_tumor_curves.pdf'],'-dpdf','-r0')
    end
    
end


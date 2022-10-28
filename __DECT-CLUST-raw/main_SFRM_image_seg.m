%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning spatial mixture of functional regression models for spectral image clustering
%
% FC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          %
% choose a regression type %
%                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%model = "PRM"; % Polynomial regression mixture
%model = "SRM"; % Spline regression mixture
model = "bSRM";% B-Spline regression mixture

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose an approach %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strategy = "FunSRM";
%strategy = "ExtraFeat-GMM";
%strategy = "kmeans";
%
K = 20;%20 % number of clusters in the data

nbr_EM_runs = 1; % algo setting

spatial_smoothing = 0;
%% data (spectral image)

%for patient_name = "subject8_tumor"%"HNSCC10"
for patient_name = "HNSCC9"%"HNSCC10"
    % for patient_name = ["subject8_tumor","HNSCC2","HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10","HNSCC11","HNSCC12","HNSCC13","HNSCC15","HNSCC17","HNSCC18","HNSCC26"]
    
    % load(['../../data/',char(patient_name),'_GT.mat']);
    % load(['../../data/',char(patient_name),'.mat']);
    load(['data/',char(patient_name),'_GT.mat']);
    load(['data/',char(patient_name),'.mat']);
    [~,~,slic_min] = ind2sub(size(segm_vol_full),find(segm_vol_full,1,'first')); % lower slice containing a tumor
    slic_inds = [1,2,3,4,6,8];%,5,7,8];%1,2,3,4,5,6];%,
    %fn_save_pdf = '';  % if don't want to save fig as pdf
    fn_save_pdf = fullfile('results', char(patient_name));  % to save fig as pdf
    
    fig_slic = figure('units','normalized','outerposition',[0.4 0.2 .5 .7]); % Init result plot
    len_sp = ceil(length(slic_inds)/2);  %length(slic_inds)
    sp = 1;
    
    
    for slic_idx = 1:length(slic_inds)
        
        
        % Prepare subject slice
        subj_slic = subject{1}(:,:,slic_inds(slic_idx));
        gr_truth = segm_vol_full(:,:,slic_inds(slic_idx)+slic_min-1);
        
        % Select a ROI on one 2D slice of the 3D volume
        %     [col_obj, row_obj] = find(imdilate(gr_truth,strel('disk',30)));  % ROI around tumor
        %  lin_obj = find(ones(size(gr_truth)));
        lin_obj = find(imdilate(gr_truth,strel('disk',60)));  % ROI around tumor;
        [col_obj, row_obj] = ind2sub(size(gr_truth),lin_obj);
        rmin = max(min(row_obj),1); cmin = max(min(col_obj),1);
        rmax = min(max(row_obj),size(subj_slic,1)); cmax = min(max(col_obj),size(subj_slic,2));
        
        
        % Plot original image
        subplot(len_sp,4,sp);
        tumor_contour = bwboundaries(gr_truth);
        imshow(subj_slic(rmin:rmax,cmin:cmax)',[]), hold on, title("slice " + num2str(slic_inds(slic_idx)+slic_min-2)) % index showed on 3DSlicer
        plot_tumor_contour(tumor_contour, [rmin, cmin], [0,0.5,1]);
        sp = sp+1;
        
        
        % Store at line 'i' the corresponding coordinates and curves
        coord = [col_obj, row_obj];
        decay_curves = zeros(length(row_obj),21);
        for r=1:length(col_obj)
            for kev=1:21
                decay_curves(r,kev) = subject{kev}(row_obj(r),col_obj(r),slic_inds(slic_idx));
            end
        end
        %
        Y = decay_curves;
        
        %Y = zscore(Y);
        %Y = Y - ones(length(Y),1)*mean(Y,1);
        
        V =  coord; % scale coordinates to keep them in [0,1]
        %T = linspace(40,140,21);
        T = linspace(0, 1, 21);
        
        %% Uncomment to apply to Other (non-spatial) curve data sets (for algo testing)
        % dataname = 'waveform'; load(dataname); Y = waveform;
        % dataname = 'satellite';load(dataname); Y = satellite;
        % dataname = 'yeast_cellcycle';load(dataname); Y = yeast_cellcycle;
        %  dataname = 'phonemes'; load(dataname); Y = phonemes;
        % [n, m] = size(Y); T = linspace(0,1, m); V = ones(length(Y), 2);%not spatial data
        %%
        
        [n, m] = size(Y);
        
        %% data matrices
        % Spatial coordinates
        Curves.spatialcoord = V;% data.VoxelCoordinates = V;
        % Curves
        Curves.abscissas = T;% data.WavelengthLevels = T;
        Curves.ordinates =  Y;% data.ReflectanceValues = Y;
        
        for model=model
            %     for model=["PRM","SRM","bSRM"]
            
            %% SRM model specification
            
            switch(model)
                case('PRM')
                    p = 4; % polynomial regression degree
                    
                    regressionOptions.basis = 'polynomial';
                    regressionOptions.p = p;
                case 'SRM'
                    spline_order = 4; % 2: linear, 3: qudratic, 4: cubic, etc, for (b-)spline spatial regression mixture
                    nknots       = 10; % fixed number of internal knots, for (b-)spline spatial regression mixture
                    
                    regressionOptions.basis = 'spline';
                    regressionOptions.spline_order = spline_order;
                    regressionOptions.nknots = nknots;
                case('bSRM')
                    Bspline_order = 4; % 2: linear, 3: qudratic, 4: cubic, etc, for (b-)spline spatial regression mixture
                    nknots       = 10; % fixed number of internal knots, for (b-)spline spatial regression mixture
                    
                    regressionOptions.basis = 'B-spline';
                    regressionOptions.Bspline_order = Bspline_order;
                    regressionOptions.nknots = nknots;
                otherwise
                    error('unknown model type');
            end
            
            %% SRM Model fitting
            mixingOption = 'softmax';
            mixingOption = 'gaussian';
            tic;
            
            
            
            switch(strategy)
                case('FunSRM')
                    
                    [mixModel, mixStats] = learn_SRM_EM(Curves, K, mixingOption, regressionOptions, nbr_EM_runs);
                    %[mixModel, mixStats] = learn_SRMF_EM(Curves, K, mixingOption, regressionOptions, nbr_EM_runs);
                    
                    %show_SRM_results(Curves, mixModel, mixStats, model);
                    
                    % show_SRM_results_new(Curves, mixModel, mixStats, model);
                    
                    %figure, plot(T, mixStats.Muk)
                    figure, plot(mixStats.stored_loglik)
                    
                    %% MRF smoothing
                    
                    if(spatial_smoothing)
                        neighb = 1;
                        [mixStats.klas, K] = mrf_smoothing(coord, mixStats.klas, neighb);
                    end
                    
                case('ExtraFeat-GMM')
                    % construct B splines features
                    
                    [~, B] = designSRM(Curves.spatialcoord, Curves.abscissas, mixingOption, regressionOptions);
                    dimBeta = size(B,2);
                    Betas = zeros(n, dimBeta);
                    %close
                    for i=1:n
                        betai = (B'*B)\B'*Y(i,:)';
                        Betas(i,:) = betai;
                        %                                                   plot(Curves.abscissas, Curves.ordinates(i,:),'o');
                        %                                                   hold on, plot(Curves.abscissas, B*betai,'r');
                        %                                                   pause
                        %                                                   clf
                        %% MRF smoothing
                        
                        if(spatial_smoothing)
                            neighb = 1;
                            [mixStats.klas, K] = mrf_smoothing(coord, mixStats.klas, neighb);
                        end
                        
                    end
                    %
                    Betas = zscore(Betas);
                    Curves.coefficients = Betas;
                    
                    [mixModel, mixStats] = learn_SRM_EM_Gauss(Curves, K, mixingOption, regressionOptions, nbr_EM_runs);
                    
                    figure, plot(mixStats.stored_loglik)
                    %% MRF smoothing
                    
                    if(spatial_smoothing)
                        neighb = 1;
                        [mixStats.klas, K] = mrf_smoothing(coord, mixStats.klas, neighb);
                    end
                    
                    
                case('kmeans')
                    %                 [mixStats.klas, mixModel.Muk] = kmeans(Curves.ordinates, K, 'MaxIter', 500, 'Display', 'iter');
                    
                    max_iter_kmeans = 500;
                    n_tries_kmeans = 3;
                    verbose_kmeans = 1;
                    sol_km = myKmeans(Curves.ordinates, K , n_tries_kmeans, max_iter_kmeans, verbose_kmeans);
                    mixStats.klas = sol_km.klas;
                    mixModel.Muk = sol_km.muk;
                    %% MRF smoothing
                    
                    if(spatial_smoothing)
                        neighb = 1;
                        [mixStats.klas, K] = mrf_smoothing(coord, mixStats.klas, neighb);
                    end
                    
                    %[mixStats.klas, mixModel.Muk] = myKmeans(Curves.ordinates, K, 'MaxIter', 500, 'Display', 'iter');
                otherwise
                    error("unknown strategy")
            end
            
            fprintf('Elapsed time %f sec \n', toc);
            %show_SRM_results_new(Curves, mixModel, mixStats, model);
            
            %% Compute similarity scores
            
            lin_tum = find(gr_truth);  % tumor;
            klas_tum = ismember(lin_obj,lin_tum);
            dice_array = zeros(1,K);
            jacc_array = zeros(1,K);  % IoU
            
            for cl_id=1:K
                ind = mixStats.klas==cl_id;
                dice_array(cl_id) = dice(ind,klas_tum);
                jacc_array(cl_id) = jaccard(ind,klas_tum);  % = sum(ind & klas_tum) / sum(ind | klas_tum);
            end
            
            [max_sim, max_cl] = max(dice_array);
            
            
            %% show results on slices
            clr = jet(K+2);
            clr = clr(1:end-2,:);
            
            figure(fig_slic)
            subplot(len_sp,4,sp);
            imshow(subj_slic(rmin:rmax,cmin:cmax)',[])
            hold on
            %         title(model)
            title(sprintf('red: Dice = %0.3f, IoU = %0.3f', max_sim, jacc_array(max_cl)), 'FontWeight', 'normal')
            
            for cl_id=1:K
                
                ind = mixStats.klas==cl_id;
                if cl_id == max_cl
                    plot(row_obj(ind)-rmin+1,col_obj(ind)-cmin+1,'.', 'color',[1 0 0],'MarkerSize',7);
                else
                    plot(row_obj(ind)-rmin+1,col_obj(ind)-cmin+1,'.', 'color',clr(cl_id,:),'MarkerSize',7);
                end
                
            end
            plot_tumor_contour(tumor_contour, [rmin, cmin], [0.99,0.99,0.99]);
            
            sp = sp + 1;
            
        end
        
    end
    
    %%
    %     stop
    %
    %     Xt = tsne(X);
    %     figure,
    %     scatter(Xt(:,1),Xt(:,2),'bx');
    %     hold on,
    %     scatter(Xt(klas_tum==1,1),Xt(klas_tum==1,2),'ro')
    %     figure, gscatter(Xt(:,1),Xt(:,2),mixStats.klas)
    %
    
    %     % Save image in pdf
    %     if ~isempty(fn_save_pdf)
    %         set(fig_slic,'Units','Inches');
    %         pos = get(fig_slic,'Position');
    %         set(fig_slic,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    %         rdn = num2str(randi(1000));
    %         print(fig_slic,[fn_save_pdf,'_',rdn,'__',char(model),'.pdf'],'-dpdf','-r0')
    %     end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot tumor contour
function plot_tumor_contour(tumor_contour, xy_min, clr)
if ~isempty(tumor_contour)
    for c=1:length(tumor_contour)
        tum_cont = cell2mat(tumor_contour(c));
        col_tum = tum_cont(:,1)-xy_min(2)+1;
        row_tum = tum_cont(:,2)-xy_min(1)+1;
        for i=1:length(row_tum) % tumor in blue
            plot(row_tum(i),col_tum(i),'.','color',clr,'MarkerSize',10);  % light blue: [0,0.5,1], white: [0.99,0.99,0.99]
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Save more figures
% fig_slic = figure(2);
% set(fig_slic,'Units','Inches');
% pos = get(fig_slic,'Position');
% set(fig_slic,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig_slic,[fn_save_pdf,'_',rdn,'__',char(model),'_curveDetails.pdf'],'-dpdf','-r0')
%
% fig_slic = figure(3);
% set(fig_slic,'Units','Inches');
% pos = get(fig_slic,'Position');
% set(fig_slic,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig_slic,[fn_save_pdf,'_',rdn,'__',char(model),'_convergence.pdf'],'-dpdf','-r0')

% %% TO BE ADDED LATER
% model_selection = 0;
% if (~model_selection) % no model selection (fixed K)
%     [mixModel, mixStats] = learn_SRM_EM(Curves, K, mixingOption, regressionOptions, nbr_EM_runs);
% else % Select K from the data
%     current_BIC = -inf;
%     Kmin = 5; Kmax = 30;
%     Krange = Kmin:Kmax;
%     
%     BIC = zeros(1, length(Krange));
%     for K = Krange(1):Krange(end)
%         [mixModel_K, mixStats_K] = learn_SRM_EM(Curves, K, mixingOption, regressionOptions, nbr_EM_runs);
%         fprintf(1,'Number of segments K: %d | BIC %f \n', K, mixStats_K.BIC);
%         if mixStats_K.BIC>current_BIC
%             mixModel = mixModel_K;
%             mixStats = mixStats_K;
%             
%             current_BIC = mixStats_K.BIC;
%         end
%         BIC(K - Krange(1)+1) = mixStats_K.BIC;
%     end
% end
end

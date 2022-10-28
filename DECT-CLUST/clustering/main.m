%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Learning spatial mixture of functional regression models          %
%                       for spectral image clustering                     %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script implements our methods proposed in paper: 
% "Spectral image clustering on dual-energy CT scans using functional regression mixtures" 
% (link of our TechRxiv to come)
% 
% Segolene Brivet (segolene.brivet@mail.mcgill.ca)


% The general workflow is:
% 1.Here: main.m: learn model and get clustering results in 'mixstats_red' variable 
% 2. build_cluster_separ_idx: compute clustering index scores and Dice scores for several images and update 'mixstats_red' variable with these scores
% 3. build_results_table: compile results for several methods and save each of them as a .mat file
% 4. resutls_analysis: boxplot and statistical values to analyse methods

%%% Data %%%
% Data were saved as .mat file with the following specifications:
%    (code for our specific dataset is shared on 'dataset_builder' folder)
% For the ground truth:     the .mat file that is loaded in this script should contain 
%                           a 3D binary (0s and 1s) image scan section and be named 'segm_vol_full'
% For the subject scan:     the .mat file that is loaded in this script should contain 
%                           a cell array named 'subject' with 21 cells (for 21 energy levels) 
%                           and each cell contains a 3D image scan section around a tumor (aligned with ground truth scan section)
% The scan sections were build so that each slice in the 3D section contains a region of interest 
% (e.g. ground truth tumor) i.e. ground truth must have some '1s' in each slice.
%%%  %%%

clear all; 
close all; 
clc;

addpath('utils');


%% Choose specifications

%%%%%%%%%%%%%%%%%%%%%%%%%%% MACHINE and FOLDER NAMES %%%%%%%%%%%%%%%%%%%%%%
machine_type = 'local';
% machine_type = 'GPU';
% /!\ Please write in 'utils/build_patient_data' the path to data folders.

results_folder_name = 'results-SgMFR-Bspl';
% results_folder_name = 'results-SgMFR-poly';
% results_folder_name = 'results-SgMVFR-Bspl';
% results_folder_name = 'results-SsMFR';
% results_folder_name = 'results-GMM';
% results_folder_name = 'results-kmeans';

%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strategy = 'MFR'; init_kmeans_dep = 1;
% strategy = 'MVFR';
% strategy = 'kmeans';
% strategy = 'gmm';
% strategy = 'kmeans';

mixingOption = 'gaussian';
% mixingOption = 'softmax';

% model = "poly"; % Polynomial regression mixture
% model = "spl"; % Spline regression mixture
model = "Bspl";% B-Spline regression mixture

p = 3; % polynomial regression degree
spline_order = 4; % 2: linear, 3: qudratic, 4: cubic, etc, for (b-)spline spatial regression mixture
Bspline_order = 4; % 2: linear, 3: qudratic, 4: cubic, etc, for (b-)spline spatial regression mixture
nknots = 10; % fixed number of internal knots, for (b-)spline spatial regression mixture

nbr_EM_runs = 1;  % number of restart of EM

% lbda and nbK could be one value > 0, or a list of values > 0. 
% Spectro-spatial index
lbda = [0.075];  % default is 0.075
% Number of clusters
nbK = [40];      % default for MFR, MVFR, kmeans is K=40   % default for GMM is K=150

%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROI size on a slice: encode a square of dim roi_size*roi_size
roi_size = 150;    % default is roi_size=150

% Number of slices for 3D section
nb_slices = 6;   % default is 6. Can also input: nb_slices='all'; to take all available slices.

% tissue enhancement window: this re-adjusts intensity values range in the image
lvl = 150;
wdw = 700;

%%%%%%%%%%%%%%%%%%%%%%%%%%% CHOOSE PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_slic_subplot = 6;   % max nb of subplots on the same figure (1 slice per subplot)

show_slices = 'middle';  % if nb_slices > max nb subplots, which slices are we plotting?
% show_slices = 'bottom';
% show_slices = 'top';

plot_clusterCurves = 0;  % generate a subfigure for each cluster, and plot all curves belonging to this cluster with mean curve and std.

plot_initResults = 0;
plot_filter3 = 1;  % recommended, add a smoothing filter on the result
plot_filter5 = 0;
plot_kmeans = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Start of loops

load_saved_mdl = 0;
for patient_name = ["HNSCC2","HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10"]


%% %%%%%%%%%%%%%%%%% TO LOAD A SAVED MODEL %%%%%%%%%%%%%%%%%%%%%
% % if model has already been run and saved, 
% % switch to 1 the variable 'load_saved_mdl'. This will plot results while avoiding re-computing the model.
% % Please uncomment this section ending with %%, and comment the 'for patient_name = ...' line before.
% load_saved_mdl = 1;
% 
% fls = dir(fullfile(results_folder_name,'/HNSCC48_245__ROI150_cl40_lbda0.075_mixstats_red.mat')); % adapt this to desired model name.
% % fls can be a list of several model names.
% for nm = 1:length(fls)
%     mixstatspath = fullfile(results_folder_name,fls(nm).name);
%     mixmodelpath = [mixstatspath(1:end-17),'mixmodel.mat'];  % (1:end-12)
%     str = split(fls(nm).name,'_');
%     patient_name = str{1}; rdn = str{2};
%     roi_size = str{4}; roi_size = str2num(roi_size(4:end));
%     K = str{5}; nbK = str2num(K(3:end));
%     lbda = str{6}; lbda = str2num(lbda(5:end));
% 
% % % special case for loading GMM
% % fls = dir(fullfile(results_folder_name,'/HNSCC21_*__gmm_ROI150_cl150_mixstats_red.mat'));
% % for nm = 1:length(fls)
% %     mixstatspath = fullfile(results_folder_name,fls(nm).name);
% %     mixmodelpath = [mixstatspath(1:end-17),'mixmodel.mat'];
% %     str = split(fls(nm).name,'_');
% %     patient_name = str{1}; rdn = str{2};
% %     roi_size = str{5}; roi_size = str2num(roi_size(4:end));
% %     K = str{6}; K = str2num(K(3:end));

%% end of LOAD A SAVED MODEL



for lambda = lbda
    
for K = nbK
    
    
    %% Data (spectral image)
    
%     close all; 
    
    patient_nm = char(patient_name); disp(patient_nm);
%     fn_save_pdf = '';  % if don't want to save plot figures, encode empty string
    fn_save_pdf = fullfile(results_folder_name, patient_nm);  % to save fig as .fig and as .pdf
    
    
    %%% Build dataset for a specific patient
    res = build_patient_data(machine_type, patient_nm, nb_slices, roi_size, max_slic_subplot, show_slices);
    
    Y = res.Y; decay_curves = res.decay_curves; gr_truth = res.gr_truth;
    klas_tum = res.klas_tum; lin_obj = res.lin_obj; subj_slic = res.subj_slic; focus = res.focus;
    slic_show = res.slic_show; slic_min_idx = res.slic_min_idx; slic_min = res.slic_min; slic_inds = res.slic_inds;
    rmin = res.rmin; rmax = res.rmax; cmin = res.cmin; cmax = res.cmax; coord = res.coord;
    
    [n, m] = size(Y);
    V = coord; V(:,3) = coord(:,3)*2;   % spacing in Z is 1.25mm, spacing in X and in Y is 0.61mm, so Z is 2*bigger than X-Y.
    T = linspace(0, 1, 21);
    
    Curves.spatialcoord = V;
    Curves.abscissas = T;
    Curves.ordinates =  Y;
    
    %% Plot original image with tissue enhancement
    
    subj_slic(subj_slic<(lvl-(wdw/2))) = lvl-(wdw/2);
    subj_slic(subj_slic>(lvl+(wdw/2))) = lvl+(wdw/2);
    
    len_sp = min(ceil(length(slic_show)/2),ceil(max_slic_subplot/2));
    slics = subj_slic(rmin:rmax,cmin:cmax,:);
    
    for i=1:length(slic_show)
        tumor_contour_list{i} = bwboundaries(gr_truth(:,:,slic_show(i)));
    end
    
    % original clustering
    if plot_initResults
        fig_slic = figure('units','normalized','outerposition',[0 0 1 1]);
        for i=1:length(slic_show)
            if i > max_slic_subplot, break; end
            subplot(len_sp,4,2*i-1);
            imshow(slics(:,:,slic_show(i)),[]), hold on, title("slice " + num2str(slic_inds(1)+slic_show(i)-1+slic_min-2)) % index showed on 3DSlicer
            plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
        end
    end
    
    % original clustering with filter
    if plot_filter3
        fig_slic2 = figure('units','normalized','outerposition',[0 0 1 1]); 
        for i=1:length(slic_show)
            if i > max_slic_subplot, break; end
            subplot(len_sp,4,2*i-1);
            imshow(slics(:,:,slic_show(i)),[]), hold on, title("slice " + num2str(slic_inds(1)+slic_show(i)-1+slic_min-2)) % index showed on 3DSlicer
            plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
        end
    end
    
    % original clustering with other filter
    if plot_filter5
        fig_slic3 = figure('units','normalized','outerposition',[0 0 1 1]); 
        for i=1:length(slic_show)
            if i > max_slic_subplot, break; end
            subplot(len_sp,4,2*i-1);
            imshow(slics(:,:,slic_show(i)),[]), hold on, title("slice " + num2str(slic_inds(1)+slic_show(i)-1+slic_min-2)) % index showed on 3DSlicer
            plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
        end
    end
    
    if plot_kmeans
        fig_slic6 = figure('units','normalized','outerposition',[0 0 1 1]);  
        for i=1:length(slic_show)
            if i > max_slic_subplot, break; end
            subplot(len_sp,4,2*i-1);
            imshow(slics(:,:,slic_show(i)),[]), hold on, title("slice " + num2str(slic_inds(1)+slic_show(i)-1+slic_min-2)) % index showed on 3DSlicer
            plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
        end
    end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MODEL FITTING                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

    if load_saved_mdl

        switch(strategy)
            case('kmeans')
                img_4D = zeros(rmax-rmin+1, cmax-cmin+1, length(slic_inds), 3);
                kev_idx = 1;
                for kev=[1,11,21]
                    img_4D(:,:,:,kev_idx) = single(permute(map_to_0_1(res.subject{kev}(cmin:cmax,rmin:rmax,slic_inds)), [2 1 3]));  % .*focus
                    kev_idx = kev_idx+1;
                end
                load(mixstatspath);
                
            case('gmm')
                load(mixstatspath);
%                 load(mixmodelpath);
                mixstats.klas = mixstats_red.klas;
%                 mixstats.Muk = mixmodel.mu;
                
            otherwise
                % if use mixstats
                % load(mixmodelpath);
                % load(mixstatspath);
                
                % if use mixstats_red
                load(mixstatspath);
                mixstats.klas = mixstats_red.klas;
                mixstats.Muk = mixstats_red.Muk;
        end

    else

        % Model specification
        switch(model)
            case('poly')
                regressionOptions.basis = 'polynomial';
                regressionOptions.p = p;
            case('spl')
                regressionOptions.basis = 'spline';
                regressionOptions.spline_order = spline_order;
                regressionOptions.nknots = nknots;
            case('Bspl')
                regressionOptions.basis = 'B-spline';
                regressionOptions.Bspline_order = Bspline_order;
                regressionOptions.nknots = nknots;
            otherwise
                error('unknown model type');
        end
        
        %% Start learning
        
        disp('Starting EM...')
        tic;
        switch(strategy)

            case('MFR')
                
                nb_EM = 1; init_kmeans = init_kmeans_dep;
                while nb_EM > 0 && nb_EM < 7  % restart if fails in converging the first time
                    try
                        [mixmodel, mixstats] = learn_SRM_EM(Curves, K, mixingOption, regressionOptions, nbr_EM_runs, lambda, init_kmeans);
                        disp(['Optim successful after ',num2str(nb_EM),' EM run(s).']);
                        nb_EM = 0;
                    catch
                        disp("Optim failed, try again.");
                        nb_EM = nb_EM + 1;
                        init_kmeans = 0;
                    end
                end
                if nb_EM == 7
                    disp("Switch to next patient.");
                    continue;
                end

            case('MVFR')
                % construct regression features
                [~, B] = designSRM(Curves.spatialcoord, Curves.abscissas, mixingOption, regressionOptions);
                dimBeta = size(B,2);
                Betas = zeros(n, dimBeta);
                for i=1:n
                    betai = (B'*B)\B'*Y(i,:)';
                    Betas(i,:) = betai;
                end
                %
                Betas = zscore(Betas);
                Curves.coefficients = Betas;

                nb_EM = 1;
                while nb_EM > 0 && nb_EM < 7  % restart if fails in converging the first time
                    try
                        [mixmodel, mixstats] = learn_SRM_EM_Gauss(Curves, K, mixingOption, regressionOptions, nbr_EM_runs, lambda);
                        disp(['Optim successful after ',num2str(nb_EM),' EM run(s).']);
                        nb_EM = 0;
                    catch
                        disp("Optim failed, try again.");
                        nb_EM = nb_EM + 1;
                    end
                end
                if nb_EM == 7
                    disp("Switch to next patient.");
                    break;
                end
                
                mixstats.Muk = mixmodel.Muk';
                

            case('kmeans')

                img_4D = zeros(rmax-rmin+1, cmax-cmin+1, length(slic_inds), 3);
                kev_idx = 1;
                for kev=[1,11,21]
                    img_4D(:,:,:,kev_idx) = single(permute(map_to_0_1(res.subject{kev}(cmin:cmax,rmin:rmax,slic_inds)), [2 1 3]));  % .*focus
                    kev_idx = kev_idx+1;
                end
                [X,Y,Z] = meshgrid(1:size(img_4D,2),1:size(img_4D,1),(1:size(img_4D,3)));
                % featureSet = cat(4,img_4D, X, Y, Z);
                featureSet = cat(4,img_4D, imgaussfilt3(single(img_4D(:,:,:,1)),1), X, Y, Z);

                mixstats_red.T = imsegkmeans3(featureSet, K, 'NormalizeInput',true);

                
            case('gmm')
                mixmodel = fitgmdist(Curves.ordinates,K,'RegularizationValue',0.1,...
                                     'Options',statset('MaxIter',1500,'TolFun',1e-6));
                [mixstats.klas, mixstats.negloglik, mixstats.posterior, mixstats.logpdf, mixstats.mahadist] = ...
                                                                                  cluster(mixmodel,Curves.ordinates);
                
                
            otherwise
                error("unknown strategy")
        end

        elaps_time_clust = toc;
        fprintf('Elapsed time for clustering algo: %f sec \n', elaps_time_clust);


    end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               RESULTS                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Colors and variables
    
    clr = jet(K+round(1/3*K));
    clr = clr(1:K,:);  % this removes red color, we will keep red for tumor clusters
    
    nb_subplot = min(max_slic_subplot,length(slic_show));
    vars.slic_show = slic_show;
    vars.slic_inds_1 = slic_inds(1);
    vars.subj_slic = subj_slic;
    vars.len_sp = len_sp;
    vars.K = K;
    vars.clr = clr;
    vars.mnx = [lvl-wdw/2, lvl+wdw/2];
    vars.max_slic_subplot = max_slic_subplot;
    vars.rmin = rmin;
    vars.rmax = rmax;
    vars.cmin = cmin;
    vars.cmax = cmax;
    vars.tumor_contour_list = tumor_contour_list;
      
        
    %%  Build label maps  &  Compute similarity scores

    switch(strategy)

        case('kmeans')

            [max_dice, max_jacc, match_klas] = compute_best_dicjac(logical(gr_truth(rmin:rmax,cmin:cmax,:)), mixstats_red.T, K);
            vars.match_klas = match_klas;
            
            if plot_kmeans
                figure(fig_slic6);
                show_result_on_img_kmeans(mixstats_red.T, img_4D, focus(rmin:rmax,cmin:cmax,slic_show), vars)
                suptitle(sprintf('%d slices / %d.   Best cluster: Dice = %0.3f, IoU = %0.3f', nb_subplot, length(slic_inds), max_dice, max_jacc))
            end
            
            if ~isempty(fn_save_pdf)
                if ~load_saved_mdl
                    rdn = num2str(randi(1000));
                    mixstats_red.dice = max_dice;
                    mixstats_red.elaps_time = elaps_time_clust;
                    save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_imsegkmeans_mixstats_red.mat'],'mixstats_red')
                end
                if plot_kmeans
                    savefig(fig_slic6,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_imsegkmeans']);
                    save_pdf(fig_slic6, [fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_imsegkmeans.pdf']);
                end
            end
            

        otherwise
        %%
        % save result variables
        
        if strcmp(strategy,'gmm')
            if ~isempty(fn_save_pdf) && ~load_saved_mdl
                rdn = num2str(randi(1000));
                save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_GMM.mat'],'mixstats')
            end
        else
            if ~isempty(fn_save_pdf)  &&  ~load_saved_mdl
                rdn = num2str(randi(1000));
                save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_lbda',num2str(lambda),'_mixmodel.mat'],'mixmodel')
                save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_lbda',num2str(lambda),'_mixstats.mat'],'mixstats')
            end
        end

        
        % Reconstruct results label map as 3D volume
        reco_lbl = zeros(size(gr_truth));
        reco_lbl(lin_obj) = mixstats.klas;
        % Simi score
        [max_dice, max_jacc, match_klas] = compute_best_dicjac(klas_tum, mixstats.klas, K);
        
        % Apply majority filter
        if plot_filter3
            fsz1 = 3; % filter size: odd nb
            reco_lbl_filt1 = medfilt3(reco_lbl, [fsz1 fsz1 fsz1]); % median filter
            % reco_lbl_filt = modefilt(reco_lbl, [fsz fsz fsz]);  % when using Matlab from 2020a
            klas_filt1 = reco_lbl_filt1(lin_obj);   % Reconstruct filtered label map as column vector
            % Simi score
            [max_dice_filt1, max_jacc_filt1, match_klas_filt1] = compute_best_dicjac(klas_tum, klas_filt1, K);
        end

        % Apply other majority filter
        if plot_filter5
            fsz2 = 5; % filter size: odd nb
            reco_lbl_filt2 = medfilt3(reco_lbl, [fsz2 fsz2 fsz2]); % median filter
            % reco_lbl_filt = modefilt(reco_lbl, [fsz2 fsz2 fsz2]);  % when using Matlab from 2020a
            klas_filt2 = reco_lbl_filt2(lin_obj);   % Reconstruct filtered label map as column vector
            % Simi score
            [max_dice_filt2, max_jacc_filt2, match_klas_filt2] = compute_best_dicjac(klas_tum, klas_filt2, K);
        end
        
        
        %% Plot results on original image
                
        if plot_initResults
            figure(fig_slic)
            vars.match_klas = match_klas;
            show_result_on_img(reco_lbl, vars);
            sgtitle(sprintf('%d slices / %d.   %d red clusters: Dice = %0.3f, IoU = %0.3f', nb_subplot, length(slic_inds), length(match_klas), max_dice, max_jacc),'FontSize',14,'FontWeight','bold');
        end
        
        % % With majority filter
        if plot_filter3
            figure(fig_slic2)
            vars.match_klas = match_klas_filt1;
            show_result_on_img(reco_lbl_filt1, vars);
            sgtitle(sprintf('%d slices / %d.   With %d*%d*%d filter  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', nb_subplot, length(slic_inds), fsz1, fsz1, fsz1, length(match_klas), max_dice_filt1, max_jacc_filt1),'FontSize',14,'FontWeight','bold');
        end
        
        % % With majority filter 2
        if plot_filter5
            figure(fig_slic3)
            vars.match_klas = match_klas_filt2;
            show_result_on_img(reco_lbl_filt2, vars);
            sgtitle(sprintf('%d slices / %d.   With %d*%d*%d filter  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', nb_subplot, length(slic_inds), fsz2, fsz2, fsz2, length(match_klas), max_dice_filt2, max_jacc_filt2),'FontSize',14,'FontWeight','bold');
        end
        
        
        
        %% Plot decay curve results
        
        if plot_clusterCurves
        % generate a subfigure for each cluster, and plot all curves belonging to this cluster with mean curve and 2*std curve.
            
            Curves.decay_curves = decay_curves;
            [fig_normCurves, fig_decayCurves, fig_loglik] = show_SRM_results_new(Curves, mixmodel, mixstats, model, match_klas, clr);
            
            if ~isempty(fn_save_pdf)
                save_pdf(fig_normCurves,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_mdlCurveDetails.pdf'])
                save_pdf(fig_decayCurves,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_decayCurveDetails.pdf'])
                save_pdf(fig_loglik,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_convergence.pdf'])
            end
        end

        
        %% Save work
        
        if ~isempty(fn_save_pdf)   % && ~load_saved_mdl

            if strcmp(strategy,'gmm')
                mixstats_red.klas = mixstats.klas;
                if exist('elaps_time_clust','var'), mixstats_red.elaps_time = elaps_time_clust; else, mixstats_red.elaps_time = NaN; end
                mixstats_red.dice = max_dice;
                if exist('max_dice_filt1','var'), mixstats_red.dice_mrg_flt = max_dice_filt1; else, mixstats_red.dice_mrg_flt = NaN; end
                save([fn_save_pdf,'_',rdn,'__gmm_ROI',num2str(roi_size),'_cl',num2str(K),'_mixstats_red.mat'],'mixstats_red')
                save([fn_save_pdf,'_',rdn,'__gmm_ROI',num2str(roi_size),'_cl',num2str(K),'_mixstats.mat'],'mixstats')
                save([fn_save_pdf,'_',rdn,'__gmm_ROI',num2str(roi_size),'_cl',num2str(K),'_mixmodel.mat'],'mixmodel')
            
                if plot_initResults, save_pdf(fig_slic, [fn_save_pdf,'_',rdn,'__gmm_ROI',num2str(roi_size),'_cl',num2str(K),'.pdf']); end
            
            else

                % save reduced variable to be smaller in size (fast transfers from GPU to local machine)
                mixstats_red.Muk = mixstats.Muk;
                mixstats_red.klas = mixstats.klas;
                mixstats_red.lambda = lambda;
                if exist('elaps_time_clust','var'), mixstats_red.elaps_time = elaps_time_clust; else, mixstats_red.elaps_time = NaN; end
                % mixstats_red.loglik = mixstats.loglik;
                mixstats_red.dice = max_dice;
                if exist('max_dice_filt1','var'), mixstats_red.dice_mrg_flt = max_dice_filt1; else, mixstats_red.dice_mrg_flt = NaN; end
                save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_lbda',num2str(lambda),'_mixstats_red.mat'],'mixstats_red')
                
                % % Save images
                % .fig
                if plot_initResults, savefig(fig_slic,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_lbda',num2str(lambda),'.fig']); end
                if plot_filter3, savefig(fig_slic2,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_lbda',num2str(lambda),'_filter',num2str(fsz1),'.fig']); end
                if plot_filter5, savefig(fig_slic3,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_lbda',num2str(lambda),'_filter',num2str(fsz2),'.fig']); end
                % .pdf
                if plot_initResults, save_pdf(fig_slic, [fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_lbda',num2str(lambda),'.pdf']); end
                if plot_filter3, save_pdf(fig_slic2,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_lbda',num2str(lambda),'_filter',num2str(fsz1),'.pdf']); end
                if plot_filter5, save_pdf(fig_slic3,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_size),'_cl',num2str(K),'_lbda',num2str(lambda),'_filter',num2str(fsz2),'.pdf']); end
            end
        end
    
    end
    

end  % end K
    
end  % end lambda

end  % end patients



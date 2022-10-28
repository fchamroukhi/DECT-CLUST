
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               Compute clustering separation indices                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segolene Brivet (segolene.brivet@mail.mcgill.ca)


% The general workflow is:
% 1. main.m: learn model and get clustering results in 'mixstats_red' variable 
% 2.Here: build_cluster_separ_idx: compute clustering index scores and Dice scores for several images and update 'mixstats_red' variable with these scores
% 3. build_results_table: compile results for several methods and save each of them as a .mat file
% 4. resutls_analysis: boxplot and statistical values to analyse methods

addpath('utils');

clear, clc

%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

machine_type = 'local';

results_folder_name = 'results-SgMFR-Bspl';
% results_folder_name = 'results-SgMFR-poly';
% results_folder_name = 'results-SgMVFR-Bspl';
% results_folder_name = 'results-SsMFR';
% results_folder_name = 'results-GMM';
% results_folder_name = 'results-kmeans';
% results_folder_name = 'results-SelSearch';

strategy = 'ours';
% strategy = 'gmm';
% strategy = 'kmeans';
% strategy = 'selsearch';

nb_slices = 6;      % default is 6. Can also input: nb_slices='all'; to take all available slices
max_slic_subplot = 6;
show_slices = 'middle';
roi_size = 150;
lambda = 0.075;

% tissue enhancement window
lvl = 150;
wdw = 700;


% Loop only on nb_K or on patient_names.
% So input a list either on list_K, or on list_patient_names.
% Otherwise, input one value per variable, e.g.:
list_K = [40];  % [30, 35, 40, 50, 60, 70]
list_patient_names = ["HNSCC2","HNSCC5","HNSCC8","HNSCC9","HNSCC10"];


% init vars:
n = max(length(list_K),length(list_patient_names));
nb_K = zeros(1,n);
ch_spat = zeros(1,n);
kl_spat = zeros(1,n);
ch_spec = zeros(1,n);
kl_spec = zeros(1,n);

db_spatec = zeros(8,n);  % row 1 to 4: different criteria for db_spat, row 5 to 8: different criteria for db_spec.

final_resuls_table = zeros(6,n);  % by row: 1.Dice 2.DB_spat 3.DB_spec 4.DBt_spat 5. DBt_spec, 6.time




%%%%%%%%%%%%%%%%%% LOAD DATA and COMPUTE CLUSTERING INDICES %%%%%%%%%%%%%%%

crit_idx = 0;

for patient_name = list_patient_names
    
    for K = list_K
        
        crit_idx = crit_idx + 1;
        
        %% Load subject scan
        
        patient_nm = char(patient_name); disp(patient_nm);
        % fn_save_pdf = '';  % if don't want to save figs
        fn_save_pdf = fullfile(results_folder_name, patient_nm);  % to save fig as .fig and .pdf
        
        res = build_patient_data(machine_type, patient_nm, nb_slices, roi_size, max_slic_subplot, show_slices);
        
        Y = res.Y; decay_curves = res.decay_curves; gr_truth = res.gr_truth;
        klas_tum = res.klas_tum; lin_obj = res.lin_obj; subj_slic = res.subj_slic; focus = res.focus;
        slic_show = res.slic_show; slic_min_idx = res.slic_min_idx; slic_min = res.slic_min; slic_inds = res.slic_inds;
        rmin = res.rmin; rmax = res.rmax; cmin = res.cmin; cmax = res.cmax; coord = res.coord;
        
        [n, m] = size(Y);
        V = coord; V(:,3) = coord(:,3)*2;   % spacing in Z is 1.25mm, spacing in X and in Y is 0.61mm, so Z is 2*bigger than X-Y.
        T = linspace(0, 1, 21);
        
        %% To read saved result variables in a folder
        
        fls = dir(fullfile(results_folder_name,[patient_nm,'_*_mixstats_red.mat']));
        if isempty(fls)
            db_spatec(:,crit_idx) = NaN;
            final_resuls_table(:,crit_idx) = NaN;
            continue;
        end
        
        nm = 1; %for nm = 1:length(fls)
        
        mixstatspath = fullfile(results_folder_name,fls(nm).name);
        
        %% Load saved model
        
        load(mixstatspath);
        
        if strcmp(strategy,'kmeans')
            
            K = length(unique(mixstats_red.T));
            nb_K(crit_idx) = K;
            
            % % Segmentation score
            [max_dice, max_jacc, match_klas] = compute_best_dicjac(logical(gr_truth(rmin:rmax,cmin:cmax,:)), mixstats_red.T, K);
            
            % % Merge tumorous clusters
            lbl = mixstats_red.T;
            for kt=2:length(match_klas)
                it = find(mixstats_red.T == match_klas(kt));
                lbl(it) = match_klas(1);
            end
            % % remove clusters outside of focus (air data)
            lbl(~focus(rmin:rmax,cmin:cmax,:)) = 0;
            % % transpose to vector label corresponding to coordinates in V
            labels_mat = zeros(size(res.gr_truth));
            labels_mat(rmin:rmax,cmin:cmax,:) = lbl;
            elem_coord = res.coord;
            elem_coord(:,3) = elem_coord(:,3)-slic_inds(1)+1;
            labels = zeros(size(elem_coord,1),1);
            for i=1:size(elem_coord,1)
                labels(i) = labels_mat(elem_coord(i,1), elem_coord(i,2), elem_coord(i,3));
            end
            
        elseif strcmp(strategy,'selsearch')
            
            str = split(fls(nm).name,'_');
            z_slic = str{end-2}; z_slic = str2double(z_slic(end));
            
            K = length(unique(mixstats_red.T));
            nb_K(crit_idx) = K;
            
            % % Segmentation score
            [max_dice, max_jacc, match_klas] = compute_best_dicjac(logical(gr_truth(rmin:rmax,cmin:cmax,z_slic)), mixstats_red.T, K);
            
            % % Merge tumorous clusters
            lbl = mixstats_red.T;
            for kt=2:length(match_klas)
                it = find(mixstats_red.T == match_klas(kt));
                lbl(it) = match_klas(1);
            end
            % % remove clusters outside of focus (air data)
            lbl(~focus(rmin:rmax,cmin:cmax,3)) = 0;
            % % transpose to vector label corresponding to coordinates in V
            labels_mat = zeros(size(res.gr_truth));
            labels_mat(rmin:rmax,cmin:cmax,3) = lbl;
            elem_coord = res.coord;
            elem_coord(:,3) = elem_coord(:,3)-slic_inds(1)+1;
            labels = zeros(size(elem_coord,1),1);
            for i=1:size(elem_coord,1)
                labels(i) = labels_mat(elem_coord(i,1), elem_coord(i,2), elem_coord(i,3));
            end
            
            
        else  % 'ours' or 'gmm'
            if ~strcmp(strategy,'gmm')
                mean_curves = mixstats_red.Muk;
            end
            
            % % Reconstruct results label map as 3D volume
            reco_lbl = zeros(size(gr_truth));
            reco_lbl(lin_obj) = mixstats_red.klas;
            % Apply median filter
            fsz1 = 3; % filter size: odd nb
            reco_lbl_filt1 = medfilt3(reco_lbl, [fsz1 fsz1 fsz1]);
            % reco_lbl_filt = modefilt(reco_lbl, [fsz fsz fsz]);  % when using Matlab from 2020a
            klas_filt1 = reco_lbl_filt1(lin_obj);   % Reconstruct filtered label map as column vector
            
            
            K = length(unique(mixstats_red.klas));
            nb_K(crit_idx) = K;
            
            % % Segmentation score
            [max_dice, max_jacc, match_klas] = compute_best_dicjac(klas_tum, klas_filt1, K);
            
            % % Merge tumorous clusters
            labels = klas_filt1;
            for kt=2:length(match_klas)
                it = find(klas_filt1 == match_klas(kt));
                labels(it) = match_klas(1);
            end
            
        end
        
        final_resuls_table(1,crit_idx) = max_dice;
        
        %% Compute validity indices
        
        % % Spatial - general
        data = V;
        [indx,~,~,~] = valid_clusterIndex(data,labels,[]);
        
        final_resuls_table(2,crit_idx) = indx(2);   % DB_spat
        % other clustering indices: CH, KL
        ch_spat(crit_idx) = indx(3);
        kl_spat(crit_idx) = indx(4);
        
        % % Spatial - focused on tumor cluster
        db_spat_list = valid_DB_index(V,labels, match_klas(1));
        
        final_resuls_table(4,crit_idx) = nanmax(db_spat_list);  % DBt_spat
        % other criteria
        db_spatec(1,crit_idx) = nanmedian(db_spat_list);
        db_spatec(2,crit_idx) = nanmean(db_spat_list);
        db_spatec(3,crit_idx) = nanmax(db_spat_list);
        db_spatec(4,crit_idx) = nanstd(db_spat_list);
        
        
        % % Curves - general
        data = Y;
        if ~strcmp(strategy,'ours') || size(mean_curves,1) ~= 21
            [indx,ssw,sw,sb] = valid_clusterIndex(data,labels,[]);
        else
            [indx,ssw,sw,sb] = valid_clusterIndex(data,labels,mean_curves');
        end
        final_resuls_table(3,crit_idx) = indx(2);        % DBt_spat
        % other clustering indices: CH, KL
        ch_spec(crit_idx) = indx(3);
        kl_spec(crit_idx) = indx(4);
        
        % % Curves - focused on tumor cluster
        db_spec_list = valid_DB_index(Y,labels, match_klas(1));
        
        final_resuls_table(5,crit_idx) = nanmax(db_spec_list);   % DBt_spec
        % other criteria
        db_spatec(5,crit_idx) = nanmedian(db_spec_list);
        db_spatec(6,crit_idx) = nanmean(db_spec_list);
        db_spatec(7,crit_idx) = nanmax(db_spec_list);
        db_spatec(8,crit_idx) = nanstd(db_spec_list);
        
        
        % % Time
        final_resuls_table(6,crit_idx) = mixstats_red.elaps_time;
        
        
        mixstats_red.dice = final_resuls_table(1,crit_idx);
        mixstats_red.DB_spat = final_resuls_table(2,crit_idx);
        mixstats_red.DB_spec = final_resuls_table(3,crit_idx);
        mixstats_red.DBt_spat = final_resuls_table(4,crit_idx);
        mixstats_red.DBt_spec = final_resuls_table(5,crit_idx);
        
        save(mixstatspath, 'mixstats_red')
        
        
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % mean of overall results % % %

mean_dice = nanmean(final_resuls_table(1,:))
mean_DBspat = nanmean(final_resuls_table(2,:))
mean_DBspec = nanmean(final_resuls_table(3,:))
mean_DBt_spat = nanmean(final_resuls_table(4,:))
mean_DBt_spec = nanmean(final_resuls_table(5,:))
mean_time = nanmean(final_resuls_table(6,:))
mean_nb_K = nanmean(nb_K)


% % % Print results one column at a time % % %
crit_idx = 1; %5
sprintf('Dice = %0.3f \n DB = %0.2f / %0.2f \n DB_t = %0.2f / %0.2f \n time = %0.2f s', ...
    final_resuls_table(1,crit_idx), final_resuls_table(2,crit_idx), final_resuls_table(3,crit_idx), final_resuls_table(4,crit_idx), final_resuls_table(5,crit_idx), final_resuls_table(6,crit_idx))



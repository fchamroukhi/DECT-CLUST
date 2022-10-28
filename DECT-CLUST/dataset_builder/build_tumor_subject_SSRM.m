

%% LOAD DATASET AND GROUND TRUTH SEGMENTATION

% This code reads DICOM files, build a subject scan section around a tumor,
% and save it as a .mat file.
% For the ground truth: a 3D binary image is saved ('segm_vol_full' variable)
% For the subject scan: a 'subject' variable is a cell array where each
% cell specifies an energy level (21 levels between 40keV and 140keV) and 
% contains a 3D image scan section around a tumor (aligned with ground truth scan section)

% Segolene Brivet (segolene.brivet@mail.mcgill.ca)


% close all; clc
% clear all;

addpath('utils');

for patient_name = ["HNSCC9"]

% for patient_name = ["HNSCC2","HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10",...
%         "HNSCC11","HNSCC12","HNSCC13","HNSCC15","HNSCC15A","HNSCC17","HNSCC17A","HNSCC18","HNSCC20",...
%         "HNSCC21","HNSCC22","HNSCC22A","HNSCC25","HNSCC26","HNSCC27","HNSCC29","HNSCC30",...
%         "HNSCC31A","HNSCC32","HNSCC33","HNSCC34","HNSCC35","HNSCC36","HNSCC37A","HNSCC38","HNSCC39",...
%         "HNSCC41","HNSCC42","HNSCC44","HNSCC44AM","HNSCC45","HNSCC46","HNSCC47","HNSCC48","HNSCC49",...
%         "HNSCC51","HNSCC52","HNSCC52AM","HNSCC53","HNSCC55","HNSCC56","HNSCC57",...
%         "HNSCC61A","HNSCC62","HNSCC63","HNSCC63A","HNSCC64A","HNSCC65A","HNSCC66","HNSCC67","HNSCC68","HNSCC69",...
%         "HNSCC70A","HNSCC71","HNSCC72A","HNSCC73","HNSCC74","HNSCC75","HNSCC76","HNSCC77","HNSCC78","HNSCC79","HNSCC80",...
%         "HNSCC81","HNSCC82","HNSCC83","HNSCC84","HNSCC85","HNSCC87","HNSCC88","HNSCC89","HNSCC90",...
%         "HNSCC91","HNSCC92","HNSCC95","HNSCC96","HNSCC97","HNSCC98",...
%         "HNSCC100","HNSCC101","HNSCC103","HNSCC105","HNSCC106","HNSCC108","HNSCC109"]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data_folder = "C:\\Users\\Segolene\\Documents\\Canada\\McGill\\PhD\\Multi-energy CT\\data\\HNSCC";
% segm_folder = "C:\\Users\\Segolene\\Documents\\Canada\\McGill\\PhD\\Multi-energy CT\\data\\HNSCC";
data_folder = "~/Documents/Data/";
segm_folder = "~/Documents/Data/";
% data_folder = "/Users/Shared/datasts/HNSCC/Multi-energy";
% segm_folder = "/Users/Shared/datasts/HNSCC/Multi-energy";
segm_type = 'tumor';
kev_list = 40:5:140;
verbose = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[char_subj, subject_nb, subject_id, ~] = get_subject_names(patient_name, segm_type, {});
study_type = '';
medic_img_struct = '';


%%% LOAD DATASET

% Load GT contour
gt_path_name = fullfile('..','data_tumor',[subject_id,'_GT.mat']);
if isfile(gt_path_name)
    segm_vol_full = load(gt_path_name).segm_vol_full;
    if verbose
        disp("Ground truth '" + gt_path_name + "' has been loaded."); disp(' ');
    end
else
    if strcmp(char_subj(1:5),'HNSCC')
        [segm_vol_full, medic_img_struct] = load_hnscc_segm(segm_folder, patient_name, segm_type);
    else
        [medic_img_struct, segm_vol_full, ~, study_type] = ...
            load_organ_segm(segm_folder, patient_name, segm_type, {});
    end
    % save for next time
    save(gt_path_name,'segm_vol_full');
    if verbose
        disp("Ground truth saved with the name '" + gt_path_name + "'."); disp(' ');
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load subject image
subject_path_name = fullfile('..','data_tumor',[subject_id,'.mat']);

if isfile(subject_path_name)
    subject = load(subject_path_name).subject;
    if verbose
        disp("Patient " + patient_name + " loaded."); disp(' ');
    end
    
    % find slices to work on
    lin_ind = find(segm_vol_full);
    [~,~,z] = ind2sub(size(segm_vol_full),lin_ind);
    slices_range = (min(z):max(z))'; % in Matlab indices (then subtract 1 for Slicer indices)
    % /!\ if ever we compute segm_vol_full and the slices are misaligned
    % with the subject section loaded in here, there would be no way to 
    % figure this out and the slices_range will be misaligned with subject
    % /!\
else
    
    if strcmp(char_subj(1:5),'HNSCC')
        subject_full = build_patient_img_hnscc(data_folder, patient_name, kev_list);
    else
%         % in case GT was saved but not subject... the study type would be unkown
%         if ~exist('study_type','var')     % 'study_type' var exists already in this script
%             [~, ~, ~, study_type] = load_organ_segm(segm_folder, patient_name, segm_type, {});
%         end
        subject_full = build_patient_img_srg(data_folder, patient_name, study_type);
    end
    
    % if GT and subject don't have the same nb of slices,
    % either duplicate some GT slices to fit subject slices, or select some slices in GT to match the fewer in subject
    % /!\ This is not the ideal solution /!\
    if size(segm_vol_full,3) ~= size(subject_full{1},3)
        
        [select_ind,vrb] = adjust_GT_subj_slices(size(segm_vol_full,3), size(subject_full{1},3));
        
        switch(vrb)
            case 'subj'
                for kev=1:length(subject_full)
                    subject_full{kev} = subject_full{kev}(:,:,select_ind);   % HNSCC60, max select_ind = 259, must not exceed 256.
                end
                    
            case 'gt'      
                for c=1:length(segm_vol_full_list)
                    segm_vol_full = segm_vol_full(:,:,select_ind);
                    
                    % rewrite ground truth
                    save(gt_path_name,'segm_vol_full');
                    if verbose
                        disp("Ground truth '" + gt_path_name + "' has been updated."); disp(' ');
                    end
                end
            otherwise
                % never happens
        end
        
    end

    % find slices to work on
    lin_ind = find(segm_vol_full);
    [~,~,z] = ind2sub(size(segm_vol_full),lin_ind);
    slices_range = (min(z):max(z))'; % in Matlab indices (then subtract 1 for Slicer indices)

    subject = cell(1, length(kev_list));
    for kev=1:length(kev_list)
        subject{kev} = subject_full{kev}(:,:,slices_range);  % slices_range-1 ?? nope
    end
    % save for next time
    save(subject_path_name,'subject');
    if verbose
        disp("Scan section saved with the name '" + subject_path_name + "'."); disp(' ');
    end
end

% close all;
clear; %clc;

end

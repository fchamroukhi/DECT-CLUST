
function subject = build_patient_img_srg(data_folder, patient_name, study_type)

disp("Loading patient " + patient_name + "...")

%%% LOAD DATA INFO %%%
dicomdir_path = fullfile(data_folder, patient_name, "DICOMDIR");
data_info = parseDicomdir(dicomdir_path);

% PixelSpacing, SliceThickness
% 

%%% GET ENERGY VALUES LIST (kev) %%%
% 2 records of energy are sometimes in files: 'Mono' and maybe 'OM' (Open Mouth)
[kev_list_Mono, kev_list_OM] = get_kev_lists(data_info);
if (isempty(kev_list_Mono) && strcmp(study_type,'Mono')) ...
        || (isempty(kev_list_OM) && strcmp(study_type,'OM'))
    error("Could not find the study corresponding to the segmentation.")
end


%%% READ ALL DATASET FOR SELECTED PATIENT %%%
if ~isempty(kev_list_Mono), energy_levels_count = length(kev_list_Mono);
else, energy_levels_count = length(kev_list_OM); end
subject = cell(1, energy_levels_count);

elaps_time = [];
try
    for energy_idx=1:length(data_info{1}.study{1}.series)
        tic

        % get energy keV in metadata
        descr = split(data_info{1}.study{1}.series{energy_idx}.info.SeriesDescription);
        for i=1:length(descr)
            k = str2double(descr{i}); if ~isnan(k), kev = k; end
        end

        % get path to slices folder (dicom images)
        slice_rel_path = split(data_info{1}.study{1}.series{energy_idx}.image{1}.info.ReferencedFileID,'\');
        slices_path = fullfile(data_folder, patient_name, slice_rel_path{1:end-1});

        % build 3D volume and store it in list (at increasing keV order)
        if strcmp(descr{1},'Mono') && strcmp(study_type,'Mono')
    %         subject_cut{find(kev_list_Mono==kev)} = read_DICOM_from_dir(slices_path, slices_range);
            subject{find(kev_list_Mono==kev)} = read_dicom(slices_path);
            stop = toc;
            elaps_time = [elaps_time; stop];
            fprintf("3D volume for energy %dkeV %s has been created in %.2fs.\n", kev, descr{1}, stop)
        elseif ~strcmp(descr{1},'Mono') && strcmp(study_type,'OM')
    %         subject_cut{find(kev_list_OM==kev)} = read_DICOM_from_dir(slices_path, slices_range);
            subject{find(kev_list_OM==kev)} = read_dicom(slices_path);
            stop = toc;
            elaps_time = [elaps_time; stop];
            fprintf("3D volume for energy %dkeV %s has been created in %.2fs.\n", kev, descr{1}, stop)
        end
    
    end
    disp("Done.")
    disp(" ")
    
catch
    disp(['Got ',num2str(energy_idx-1),' energy levels in subject. Failed at ',num2str(energy_idx),'. Stop loading...']);
end

end


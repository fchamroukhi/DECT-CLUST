
function subject = build_patient_img_hnscc(data_folder, patient_name, kev_list)

disp("Loading patient " + patient_name + "...")

pn = char(patient_name);
subject = cell(1, length(kev_list));
elaps_time = [];

fn_list = dir(fullfile(data_folder,patient_name,'CT_*kev*'));
if length(kev_list)~=length(fn_list)
    error('Nb of folder differs from nb of energy levels');
end
for nm = 1:length(fn_list)
    tic
    kev_folder = fn_list(nm).name;
    image_folder = fullfile(data_folder,patient_name,kev_folder, ['Hnscc',pn(6:end)]);
    kev_cell = regexp(kev_folder,'_(\d+)kev','tokens');
    kev = (str2double(kev_cell{1}{:})-35)/5;
    
% for kev=1:length(kev_list)
%     tic
%     image_folder = fullfile(data_folder,patient_name,['CT_',num2str(kev_list(kev)),'kev'], ['Hnscc',pn(6:end)]);
    subject{kev} = read_DICOM_from_dir(image_folder, []);
    stop = toc;
    elaps_time = [elaps_time; stop];
    fprintf("3D volume for energy %dkeV has been created in %.2fs.\n", kev_list(kev), stop)
end

disp("Done.")
disp(" ")

end


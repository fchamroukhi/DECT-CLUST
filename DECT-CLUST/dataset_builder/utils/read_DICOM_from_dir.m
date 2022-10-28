
function [imaVOL] = read_DICOM_from_dir(inputdir, slice_range)
% imaVOL is a 3D volume of slices
% inputdir is the folder containing all dicom slices
% slice_range is a vector of slices to be read

dname=dir(fullfile(inputdir,'IM*'));

% if ~isempty(slice_range)
%     % selected slices to read (in the order they appear)
%     sl=slice_range;
% else

% all slices to read
sl=zeros(1,length(dname));
for i = 1:length(dname)
    dcm = fullfile(inputdir, dname(i).name);
    di = dicominfo(dcm);
%         sl(i) = di.InstanceNumber;
    sl(i) = di.ImagePositionPatient(3);
end
sl=sort(sl);
    
% end

% Read slices in the order of sl
taken_sl = [];
init_mat = 0;
for i = 1:length(dname)
    dcm = fullfile(inputdir, dname(i).name);
    di = dicominfo(dcm);
%     if ismember(di.InstanceNumber, sl)
    if ismember(di.ImagePositionPatient(3), sl)
        d = int16(dicomread(dcm));
        if(init_mat==0)
            imaVOL = zeros(size(d,2),size(d,1),length(sl));
            init_mat = 1;
        end
        dr = d*di.RescaleSlope + di.RescaleIntercept;

        dr = double(dr);
%         imaVOL(:,:,find(sl==di.InstanceNumber)) = dr';
%         taken_sl = [taken_sl; di.InstanceNumber];
        imaVOL(:,:,find(sl==di.ImagePositionPatient(3))) = dr';
        taken_sl = [taken_sl; di.ImagePositionPatient(3)];
    end
end

missing = setdiff(sl,taken_sl);
if length(missing)==length(sl)
    disp("No slice has been found nor read.")
elseif ~isempty(missing)
    fprintf("These slices have not been found: ["); fprintf('%g ', missing); fprintf(']\n');
end

end

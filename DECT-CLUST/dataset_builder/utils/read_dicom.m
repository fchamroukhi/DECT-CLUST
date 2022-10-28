
function imaVOL = read_dicom(inputdir)
% inputdir  : path of the slices for one energy level

dname=dir(fullfile(inputdir,'IM*'));

ipp3=zeros(1,length(dname));

for i=1:length(dname)
    dcm = fullfile(dname(i).folder,dname(i).name);
    d=int16(dicomread(dcm));
    di=dicominfo(dcm);
    if(i==1)
        imaVOL=zeros(size(d,2),size(d,1),length(dname));
    end
    dr=d*di.RescaleSlope + di.RescaleIntercept;
    imaVOL(:,:,i) = dr';
    ipp3(i) = di.ImagePositionPatient(3);
end

[~,index]=sort(ipp3);
imaVOL=imaVOL(:,:,index);

end

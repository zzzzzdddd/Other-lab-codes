
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
addpath('T:\Imaging\Multimodal\MRF\Peter\Joon_code\matlab_joon');

for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    head = load_nii('MNI_ROI_final.nii');
    head = head.hdr;
    files = ["sMNI_CSF.nii" "sMNI_extension.nii" "sMNI_GM" "sMNI_junction" "sMNI_thickness" "sMNI_WM"];
    for f = files
        file = load_nii(char(f));
        file.hdr = head;
        save_nii(file,char(f));
    end
end
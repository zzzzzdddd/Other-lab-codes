% subjID = [XXXX];
cd T:\Imaging\Multimodal\MRF\Peter
load('MRFvalues_withcsf.mat');
subjID = [XXXXX];

MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path);
    copyfile MNI_GM_prob.nii MNI_MRF_junc.nii;

    a = load_untouch_nii('MNI_GM_prob.nii');
    b = load_untouch_nii('MNI_WM_prob.nii');
    c = load_untouch_nii('MNI_MRF_junc.nii');

   
end
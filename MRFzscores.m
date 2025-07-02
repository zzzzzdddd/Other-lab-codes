cd T:\Imaging\Multimodal\MRF\Peter
load('MRFvalues_withcsf.mat');
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
subjID = [XXXX];
for p = subjID
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    copyfile MNI_T1.nii MNI_T1z.nii
    copyfile MNI_T2.nii MNI_T2z.nii
    t1z = load_untouch_nii('MNI_T1z.nii');
    t1zi = single(t1z.img);
    t2z = load_untouch_nii('MNI_T2z.nii');
    t2zi = single(t2z.img);
    t1z.img = (t1zi-mT1)./sdT1;
    t2z.img = (t2zi-mT2)./sdT2;
    t1z.img(isinf(t1z.img)|isnan(t1z.img)) = 0;
    t2z.img(isinf(t2z.img)|isnan(t2z.img)) = 0;
    save_untouch_nii(t1z, 'MNI_T1z.nii');
    save_untouch_nii(t2z, 'MNI_T2z.nii');
end


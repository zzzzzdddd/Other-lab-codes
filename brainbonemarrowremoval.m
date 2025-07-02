clear all
% subjID=[XXX]



MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
addpath('T:\Imaging\Multimodal\MRF\Peter\Joon_code\matlab_joon');

%% 
for p = subjID
    path = strcat(MRF_path,'\',p,'\reorient_all');
    cd(path)
    GM = load_untouch_nii('c1T1.nii');
    WM = load_untouch_nii('c2T1.nii');
    
    mask = load_untouch_nii('Brain_Mask.nii');
    mask.img = GM.img + WM.img;
%     mask.img(find(mask.img>=0.96)) = 1;
%     mask.img(find(mask.img<0.96)) = 0;
    save_untouch_nii(mask,'Brain_Mask.nii')
end

%% niftiread and write
for p = subjID
    path = strcat(MRF_path,'\',p,'\reorient_all\fatT1');
    cd(path)
%     GM = niftiread('c1T1.nii');
%     WM = niftiread('c2T1.nii');
%     copyfile T1.nii Brain_Mask.nii;
    mask = load_untouch_nii('T1_brain_mask.nii');
    image = load_untouch_nii('T1.nii');
    maski = single(mask.img);
    imagei = single(image.img);
    imagei = imagei.*maski;
    image.img = imagei;
    save_untouch_nii(image,'T1.nii');
end
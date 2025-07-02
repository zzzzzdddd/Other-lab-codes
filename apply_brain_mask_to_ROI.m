
clear all
subjID=[XXX];


MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
addpath('T:\Imaging\Multimodal\MRF\Peter\Joon_code\matlab_joon');

for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    copyfile MNI_ROI.nii MNI_ROI_final.nii;
    data = load_untouch_nii('MNI_ROI_final.nii');
        
    cd(strcat(MRF_path,'\',p,'\MRF_VBM'));
    mask = load_untouch_nii('MNI_Brain_Mask.nii');
    newROI = data.img .* mask.img;
    
    cd(path);
    im = single(newROI);
    [mxv,idx] = max(im(:));
    level = graythresh(im);
    bina = zeros(size(im));
    bina(find(im>(level*mxv))) = 1;
    bina = single(bwmorph3(bina,'clean'));
    SE = strel('disk',3);
    bina = imclose(bina,SE);
    bina = imopen(bina,SE);
    
    data.img(find(bina>0)) = 1;
    data.img(find(bina==0)) = 0;
    save_untouch_nii(data,'MNI_ROI_final.nii');

end

%% dilated ROI

% subjID=[XXXX];

MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
addpath('T:\Imaging\Multimodal\MRF\Peter\Joon_code\matlab_joon');

for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    copyfile MNI_ROI.nii MNI_ROI_final.nii;
    data = load_untouch_nii('MNI_ROI_final.nii');
        
    cd(strcat(MRF_path,'\',p,'\reorient_all'));
    mask = load_untouch_nii('MNI_Brain_Mask.nii');
    newROI = data.img .* mask.img;
    
    cd(path);
    im = single(newROI);
    [mxv,idx] = max(im(:));
    level = graythresh(im);
    bina = zeros(size(im));
    bina(find(im>(level*mxv))) = 1;
    bina = single(bwmorph3(bina,'clean'));
    SE = strel('disk',3);
    bina = imclose(bina,SE);
    bina = imopen(bina,SE);
    bina = imdilate(bina,SE);
    data.img = data.img .* mask.img;

    data.img(find(bina>0)) = 1;
    data.img(find(bina==0)) = 0;
    save_untouch_nii(data,'MNI_ROI_final.nii');

end

%% make empty label files for normal controls
clear all

subjID=[XXX];

MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
addpath('T:\Imaging\Multimodal\MRF\Peter\Joon_code\matlab_joon');

for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    copyfile MNI_T1.nii MNI_ROI_final.nii;
    data = load_untouch_nii('MNI_ROI_final.nii');
        
    data.img = single(data.img).*0;
    save_untouch_nii(data,'MNI_ROI_final.nii');

end


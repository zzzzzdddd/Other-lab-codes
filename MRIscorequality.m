clear all
subjID=[XXX];

MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
% MRF_path= '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Patients';

MRF_path= '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Patients';
subjID=[XXXX];

MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';

cvGMall       = [];
cvWMall       = [];
cvcGMWMall    = [];
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    GM = load_untouch_nii('MNI_GM.nii');
    GMmask = single(GM.img);
    GMmask(GMmask>0.5) = 1;
    GMmask(GMmask<0.5) = 0;

    WM = load_untouch_nii('MNI_WM.nii');
    WMmask = single(WM.img);
    WMmask(WMmask>0.5) = 1;
    WMmask(WMmask<0.5) = 0;
    
    MRI = load_untouch_nii('MNI_MPRAGE.nii');
    image = single(MRI.img);

    intensitiesGM = image.*GMmask;
    intensitiesWM = image.*WMmask;
    meanGM = mean(intensitiesGM(intensitiesGM>0));
    meanWM = mean(intensitiesWM(intensitiesWM>0));
    stdGM = std(intensitiesGM(intensitiesGM>0));
    stdWM = std(intensitiesWM(intensitiesWM>0));
    cvGM       = repmat(stdGM./meanGM, [size(image(:)),1]); % coefficient of variation for GM signal intensities
    cvWM       = repmat(stdWM./meanWM, [size(image(:)),1]); % coefficient of variation for WM signal intensities
    cvcGMWM    = repmat(abs(meanGM-meanWM)./mean([stdGM stdWM]), [size(image(:)),1]); % standardized contrast of GM to WM intensities
    
    copyfile MNI_GM.nii MNI_cvGM.nii
    copyfile MNI_GM.nii MNI_cvWM.nii
    copyfile MNI_GM.nii MNI_cvcGMWM.nii
    
    file = load_untouch_nii('MNI_cvGM.nii');
    file.img = cvGM;
    save_untouch_nii(file, 'MNI_cvGM.nii');
    file = load_untouch_nii('MNI_cvWM.nii');
    file.img = cvWM;
    save_untouch_nii(file, 'MNI_cvWM.nii');
    file = load_untouch_nii('MNI_cvcGMWM.nii');
    file.img = cvcGMWM;
    save_untouch_nii(file, 'MNI_cvcGMWM.nii');
   
    cvGMall       = vertcat(cvGMall,stdGM./meanGM); % coefficient of variation for GM signal intensities
    cvWMall       = vertcat(cvWMall, stdWM./meanWM); % coefficient of variation for WM signal intensities
    cvcGMWMall    = vertcat(cvcGMWMall, abs(meanGM-meanWM)./mean([stdGM stdWM])); % standardized contrast of GM to WM intensities

end
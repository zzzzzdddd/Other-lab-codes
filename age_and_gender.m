%% patient predict and score
subjID=[XXX];

MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
cd('Z:\Imaging\Multimodal\MRF\Peter')
clinifo = readmatrix('Ptinfo.xlsx', 'Range','1:44');
age = clinifo(2:length(subjID)+1,1);
sex = clinifo(2:length(subjID)+1,2)/2;

for p = 1:length(subjID)
    p
    path = strcat(MRF_path,'\',subjID(p),'\MRF_VBM');
    cd(path);
    copyfile MNI_T1.nii MNI_age.nii
    copyfile MNI_T1.nii MNI_gender.nii

    imt1 = load_untouch_nii('MNI_age.nii');
    imt1.img(imt1.img>0) = age(p);
    save_untouch_nii(imt1, 'MNI_age.nii');

    imt2 = load_untouch_nii('MNI_gender.nii');
    imt2.img(imt2.img>0) = sex(p);
    save_untouch_nii(imt2, 'MNI_gender.nii');

end
%% volunteer
subjID=[XXXX];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
cd('Z:\Imaging\Multimodal\MRF\Peter')
clinifo = readmatrix('HCinfo.csv', 'Range','1:68');
age = clinifo(2:length(subjID)+1,1);
sex = clinifo(2:length(subjID)+1,3)/2;

for p = 1:length(subjID)
    p
    path = strcat(MRF_path,'\',subjID(p),'\MRF_VBM');
    cd(path);
    copyfile MNI_T1.nii MNI_age.nii
    copyfile MNI_T1.nii MNI_gender.nii

    imt1 = load_untouch_nii('MNI_age.nii');
    imt1.img(imt1.img>0) = age(p);
    save_untouch_nii(imt1, 'MNI_age.nii');

    imt2 = load_untouch_nii('MNI_gender.nii');
    imt2.img(imt2.img>0) = sex(p);
    save_untouch_nii(imt2, 'MNI_gender.nii');

end
%
clear all
subjID=[];

MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
MRF_path= '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Patients';

fileNames = {'MNI_MPRAGE.nii','MNI_junction.nii', 'MNI_extension.nii', 'MNI_thickness.nii',  'sMNI_extension.nii', 'sMNI_junction.nii', 'sMNI_thickness.nii', ... 
            'MNI_WM.nii', 'MNI_GM.nii', 'sMNI_GM.nii', 'sMNI_WM.nii', 'MNI_T1.nii', 'MNI_T2.nii', 'MNI_GM_prob.nii', 'MNI_WM_prob.nii'};
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    mask = load_untouch_nii('MNI_Brain_Mask.nii');
    for k = 1:numel(fileNames)
        fname = fileNames{k};
        in = load_untouch_nii(fname);
        in = double(in.img.*mask.img);
        cmean = mean(in(mask.img==1), [1 2 3]);
        cstd = std(in(mask.img==1), 0, [1 2 3]);
        out = (in - cmean)./cstd;
        out = out.*mask.img;
        output = make_nii(out);
        save_nii(output, strcat('zmean_',fname));
    end
end
%% Normal

subjID=[XXXX];

MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
MRF_path= '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Normal';


fileNames = {'MNI_MPRAGE.nii','MNI_junction.nii', 'MNI_extension.nii', 'MNI_thickness.nii',  'sMNI_extension.nii', 'sMNI_junction.nii', 'sMNI_thickness.nii', ... 
            'MNI_WM.nii', 'MNI_GM.nii', 'sMNI_GM.nii', 'sMNI_WM.nii', 'MNI_T1.nii', 'MNI_T2.nii', 'MNI_GM_prob.nii', 'MNI_WM_prob.nii'};
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    mask = load_untouch_nii('MNI_Brain_Mask.nii');
    for k = 1:numel(fileNames)
        fname = fileNames{k};
        in = load_untouch_nii(fname);
        in = double(in.img.*mask.img);
        cmean = mean(in(mask.img==1), [1 2 3]);
        cstd = std(in(mask.img==1), 0, [1 2 3]);
        out = (in - cmean)./cstd;
        out = out.*mask.img;
        output = make_nii(out);
        save_nii(output, strcat('zmean_',fname));
    end
end
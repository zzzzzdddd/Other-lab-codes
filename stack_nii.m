%% stack niis into one nii and move to folder
clear all
subjID=[XXX];
% dest = '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Patients/MRFupload_to_hpc_new';
% MRF_path= '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Patients';
% dest='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRFupload_wholetz';
dest='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRFupload_linreg';
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

% fileNames = {'MNI_MPRAGE.nii', 'MNI_Brain_Mask.nii', 'MNI_junction.nii', 'MNI_extension.nii', 'MNI_thickness.nii',  'sMNI_extension.nii', 'sMNI_junction.nii', 'sMNI_thickness.nii', ... 
%             'MNI_WM.nii', 'MNI_GM.nii', 'MNI_CSF.nii', 'sMNI_GM.nii', 'sMNI_WM.nii', 'sMNI_CSF.nii', 'MNI_T1.nii', 'MNI_T2.nii', 'MNI_GM_prob.nii', 'MNI_WM_prob.nii'};
fileNames = {'MNI_MPRAGE_fn.nii', 'MNI_MPRAGE_fn.nii', 'MNI_junction_fn.nii', 'MNI_extension_fn.nii', 'MNI_thickness_fn.nii', ... 
            'MNI_GM_fn.nii', 'MNI_WM_fn.nii', 'MNI_CSF_fn.nii', 'MNI_T1z.nii', 'MNI_T2z.nii', ...
            'sMNI_junction_fn.nii', 'sMNI_extension_fn.nii', 'sMNI_thickness_fn.nii', ... 
            'sMNI_GM_fn.nii', 'sMNI_WM_fn.nii', 'sMNI_CSF_fn.nii', 'sMNI_T1z.nii', 'sMNI_T2z.nii', ...
            'MNI_ROI_final.nii'};
fileNames = {'MNI_T1w.nii', 'MNI_T1.nii', 'MNI_T2.nii', 'MNI_GM_prob.nii', 'MNI_WM_prob.nii', ... 
            'MNI_CSF_prob.nii', 'MNI_MRFjuncz.nii', 'MNI_MRFT1juncz.nii', 'MNI_MRFT2juncz.nii', 'MNI_MRFextenz.nii', ...
            'sMNI_T1_nocsf.nii', 'sMNI_T2_nocsf.nii', 'MNI_GMprobz.nii', 'MNI_WMprobz.nii', 'MNI_CSFprobz.nii','MNI_ROI_final.nii'};
fileNames = {'MNI_T1w.nii', 'MNI_T1.nii', 'MNI_T2.nii', 'MNI_GM_prob.nii', 'MNI_WM_prob.nii', ... 
            'MNI_CSF_prob.nii', 'MNI_MRFjuncz.nii', 'MNI_T1z.nii', 'MNI_T2z.nii', 'MNI_MRFextenz.nii', ...
            'MNI_GMprobz.nii', 'MNI_WMprobz.nii', 'MNI_CSFprobz.nii','MNI_ROI_final.nii'};
fileNames = {'MNI_T1w.nii', 'MNI_T1.nii', 'MNI_T2.nii', 'MNI_GM_prob.nii', 'MNI_WM_prob.nii', ... 
            'MNI_CSF_prob.nii', 'MNI_MRFjuncz.nii', 'MNI_T1LMz.nii', 'MNI_T1LMz.nii', 'MNI_MRFextenz.nii', ...
            'MNI_GMprobz.nii', 'MNI_WMprobz.nii', 'MNI_CSFprobz.nii','MNI_ROI_final.nii'};
for p = 1:size(subjID,2)
    p
    path = strcat(MRF_path,'\',subjID(p),'\MRF_VBM');
%     path = strcat(MRF_path,'/',subjID(p),'/MRF_VBM');
    cd(path)
    for k = 1:numel(fileNames)
        fname = fileNames{k};
        z = load_untouch_nii(fname); 
        y(:,:,:,k) = single(z.img);
%         if k == 2
%             z = single(load_untouch_nii('MNI_GM_fn.nii').img) + single(load_untouch_nii('MNI_WM_fn.nii').img);
%             y(:,:,:,k) = single(z);
%         else    
%             z = load_untouch_nii(fname);
%             y(:,:,:,k) = single(z.img);
%         end
    end
    output = make_nii(y);
    save_nii(output, char(strcat(subjID(p),'.nii')))
    sourceFile = strcat(subjID(p),'.nii');
%     destFile = fullfile(dest, sourceFile); 
    gzip(sourceFile, dest);
end


%% for patients move image niis into folder Irene hp

subjID=[XXX];

subjID=["XXX"];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
dest = 'C:\Users\Irene\Downloads\MRF_VBM\imagesTr';

% dest = '/Users/irene/Downloads/MRF_VBM/imagesTr';
for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    sourceFile = strcat(p,'.nii');
    destFile = fullfile(dest, sourceFile); 
    copyfile(sourceFile, destFile);
end
%% for patients move roi niis into folder
subjID=[XXX];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
dest = 'C:\Users\dingz\Downloads\MRF_VBM\labelsTr';

% MRF_path= '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Patients';
dest = 'C:\Users\irene\Downloads\MRF_VBM\labelsTr';

for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    sourceFile = 'MNI_ROI_final.nii';
    destFile = fullfile(dest, strcat(p,'.nii')); 
    copyfile(sourceFile, destFile);
end
%% for normal, stack niis   V2 and "V36_13592" missing files
clear all
subjID=[XXX];

% MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
MRF_path= '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Normal';
dest = '/Volumes/eegrvw/Imaging/Multimodal/MRF/Peter/MRF_hc_tohpc'
% fileNames = {'MNI_MPRAGE.nii', 'MNI_Brain_Mask.nii', 'MNI_junction.nii', 'MNI_extension.nii', 'MNI_thickness.nii',  'sMNI_extension.nii', 'sMNI_junction.nii', 'sMNI_thickness.nii', ... 
%             'MNI_WM.nii', 'MNI_GM.nii', 'MNI_CSF.nii', 'sMNI_GM.nii', 'sMNI_WM.nii', 'sMNI_CSF.nii', 'MNI_T1.nii', 'MNI_T2.nii', 'MNI_GM_prob.nii', 'MNI_WM_prob.nii'};
fileNames = {'MNI_MPRAGE.nii', 'MNI_Brain_Mask.nii', 'MNI_junction.nii', 'MNI_extension.nii', 'MNI_thickness.nii', ... 
            'MNI_WM.nii', 'MNI_GM.nii', 'MNI_CSF.nii', 'MNI_T1.nii', 'MNI_T2.nii', 'MNI_ROI_final.nii'};
fileNames = {'MNI_T1w.nii', 'MNI_T1.nii', 'MNI_T2.nii', 'MNI_GM_prob.nii', 'MNI_WM_prob.nii', ... 
            'MNI_CSF_prob.nii', 'MNI_ROI_final.nii'};
for p = 1:size(subjID,2)
    p
%     path = strcat(MRF_path,'\',subjID(p),'\MRF_VBM');
    path = strcat(MRF_path,'/',subjID(p),'/MRF_VBM');
    cd(path)
    for k = 1:numel(fileNames)
        fname = fileNames{k};
        z(k) = load_untouch_nii(fname);
        y(:,:,:,k) = z(k).img;
    end
    output = make_nii(y);
    save_nii(output, char(strcat(subjID(p),'.nii')))
    sourceFile = strcat(subjID(p),'.nii');
    destFile = fullfile(dest, sourceFile); 
    movefile(sourceFile, destFile);
end


%% no stacking
clear all
subjID=[XXX];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
dest = 'C:\Users\eegrvw\Downloads\MRF_VBM\Training';
% MRF_path= '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Patients';

% fileNames = {'MNI_MPRAGE.nii', 'MNI_Brain_Mask.nii', 'MNI_junction.nii', 'MNI_extension.nii', 'MNI_thickness.nii',  'sMNI_extension.nii', 'sMNI_junction.nii', 'sMNI_thickness.nii', ... 
%             'MNI_WM.nii', 'MNI_GM.nii', 'MNI_CSF.nii', 'sMNI_GM.nii', 'sMNI_WM.nii', 'sMNI_CSF.nii', 'MNI_T1.nii', 'MNI_T2.nii', 'MNI_GM_prob.nii', 'MNI_WM_prob.nii'};
fileNames = {'MNI_MPRAGE.nii', 'MNI_Brain_Mask.nii', 'MNI_junction.nii', 'MNI_extension.nii', 'MNI_thickness.nii', ... 
            'MNI_WM.nii', 'MNI_GM.nii', 'MNI_CSF.nii', 'MNI_T1.nii', 'MNI_T2.nii'};

for p = 1:size(subjID,1)
    p
    path = strcat(MRF_path,'\',subjID(p),'\MRF_VBM');
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
%     mask = load_untouch_nii('MNI_Brain_Mask_binary.nii')
    for k = 1:numel(fileNames)
        fname = fileNames{k};

    end
    output = make_nii(y);
    save_nii(output, char(strcat('P', string(p),'.nii')))

    sourceFile = strcat('P', string(p),'.nii');
    destFile = fullfile(dest, sourceFile); 
    copyfile(sourceFile, destFile);
end


file = load_untouch_nii('P98_15551.nii');
image = double(file.img);

slice = image(:,:,:,16);
slice(slice > 0) = 1;
image(:,:,:,16) = slice;

file.img = image;
save_untouch_nii(file, 'P98_15551.nii')
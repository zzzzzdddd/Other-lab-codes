subjID=[XXXXX]
MRF_path='Z:\Imaging\Multimodal\Myelin\Patients\Resection_analysis';
myelinpath = 'Z:\Imaging\Multimodal\Myelin\Patients';
for p = subjID
    cd(MRF_path)
    mkdir(p)
    path = strcat(MRF_path,'\',p);
    cd(path);
    mkdir Pre
    mkdir Post
    copyfile(strcat(myelinpath,'\',p,'\T1.nii'), strcat(MRF_path,'\',p,'\Pre\T1.nii'));
    copyfile(strcat(myelinpath,'\',p,'\postop.nii'), strcat(MRF_path,'\',p,'\Post\postop.nii'));
end
%%
list = dir('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\VXX\VBM_newANN');
filenames = string({list.name});
filenames = filenames(7:end);
subjID=[XXXX]
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
for p = subjID
    dest = strcat(MRF_path,'\',p,'\','VBM_newANN');
    path = strcat(MRF_path,'\',p);
    cd(path);
    for f = filenames
        movefile(strcat(MRF_path,'\',p,'\',f), dest);
    end
end

%%
cd T:\Imaging\Multimodal\MRF\Recon_MRF_3T
copyfile MNI152_T1_1mm_brain.nii MNI152_T1_1mm_brain_mask.nii 
mnimask = load_untouch_nii('MNI152_T1_1mm_brain_mask.nii');
image = single(mnimask.img);
image(image>0) = 1;
mnimask.img = image;
save_untouch_nii(mnimask, 'MNI152_T1_1mm_brain_mask.nii');
%% grab T2 files copy to folder for normal
% Coded by JYC; edited by TYS
addpath('T:\Imaging\Multimodal\MRF\Joon\projects\ETC_function\NIfTI_20140122')
%T1overT2 images (modified 01102020)
clear all; close all;
class = 'Normal'; %'Normal' or 'Patients'
subjID=["XXXX"];
folder = fullfile('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T', class);
out_dir = fullfile('Z:\Imaging\Multimodal\Myelin', class);

for subject_form = subjID
Seqs = dir(fullfile(folder, subject_form, 'nii'));
T1_count = 1;
strcmpt = strsplit(Seqs(3).folder, '\');
subject = strcmpt{7};
for order = 1:length(Seqs)
    if ~isempty(strfind(Seqs(order).name, 'T1_MPRAGE_AX')) & isempty(strfind(Seqs(order).name, 'COR'))
        T1_MPRAGE_AX_id(T1_count, 1) = order;
        T1_count = T1_count + 1;
    end
    if ~isempty(strfind(Seqs(order).name, 'T2_SPACE_AX')) & isempty(strfind(Seqs(order).name, 'COR'))
        T2_SPACE_AX_id = order;
    end
end
T2_file = dir(fullfile(folder, subject, 'nii', Seqs(T2_SPACE_AX_id).name, '*.nii'));
copyfile(fullfile(T2_file.folder, T2_file.name), 'T2im.nii')
movefile('T2im.nii', fullfile(out_dir, subject_form))
end

%% grab T2 files copy to folder for patients
% Coded by JYC; edited by TYS
addpath('T:\Imaging\Multimodal\MRF\Joon\projects\ETC_function\NIfTI_20140122')
%T1overT2 images (modified 01102020)
clear all; close all;
class = 'Patients'; %'Normal' or 'Patients'
subjID = ["XXXX"]; 
folder = fullfile('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T', class);
out_dir = fullfile('Z:\Imaging\Multimodal\Myelin', class);

for subject_form = subjID
Seqs = dir(fullfile(folder, subject_form, 'nii'));
T1_count = 1;
strcmpt = strsplit(Seqs(3).folder, '\');
subject = strcmpt{7};
for order = 1:length(Seqs)
    if ~isempty(strfind(Seqs(order).name, 'T1_MPRAGE_AX')) & isempty(strfind(Seqs(order).name, 'COR'))
        T1_MPRAGE_AX_id(T1_count, 1) = order;
        T1_count = T1_count + 1;
    end
    if ~isempty(strfind(Seqs(order).name, 'T2_SPACE_AX')) & isempty(strfind(Seqs(order).name, 'COR'))
        T2_SPACE_AX_id = order;
    end
end
T2_file = dir(fullfile(folder, subject, 'nii', Seqs(T2_SPACE_AX_id).name, '*.nii'));
copyfile(fullfile(T2_file.folder, T2_file.name), 'T2im.nii')
movefile('T2im.nii', fullfile(out_dir, subject_form))
end

%% patient copy files 
subjID = ["XXXX"]; 
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
myelin_path = 'Z:\Imaging\Multimodal\Myelin\Patients';
for p = subjID
    path = strcat(myelin_path,'\',p);
    cd(path);
    copyfile('MNI_T1im.nii', char(strcat(p,"T1.nii")))
    copyfile('MNI_T2im.nii', char(strcat(p,"T2.nii")))
    movefile(strcat(p,"T1.nii"),'Z:\Imaging\Multimodal\Myelin\allT1');
    movefile(strcat(p,"T2.nii"),'Z:\Imaging\Multimodal\Myelin\allT2');
    copiedmask = char(strcat(p,"mask.nii"));
    maskdest = fullfile('Z:\Imaging\Multimodal\Myelin\allmasks',copiedmask);
    copyfile('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\MNI152_T1_1mm_brain.nii',maskdest)
end

%% normal copy files 
subjID=["XXX"];
myelin_path = 'Z:\Imaging\Multimodal\Myelin\Normal';
for p = subjID
    path = strcat(myelin_path,'\',p);
    cd(path);
    copyfile('MNI_T1im.nii', char(strcat(p,"T1.nii")))
    copyfile('MNI_T2im.nii', char(strcat(p,"T2.nii")))
    movefile(strcat(p,"T1.nii"),'Z:\Imaging\Multimodal\Myelin\allT1');
    movefile(strcat(p,"T2.nii"),'Z:\Imaging\Multimodal\Myelin\allT2');
    copiedmask = char(strcat(p,"mask.nii"));
    maskdest = fullfile('Z:\Imaging\Multimodal\Myelin\allmasks',copiedmask);
    copyfile('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\MNI152_T1_1mm_brain.nii',maskdest)
end
%% 
subjID=["XXX"];
myelin_path = 'Z:\Imaging\Multimodal\Myelin\Patients';
normalizedpath1 = 'Z:\Imaging\Multimodal\Myelin\nyul_normalizedT1';
normalizedpath2 = 'Z:\Imaging\Multimodal\Myelin\nyul_normalizedT2';

for p = subjID
    path = normalizedpath1;
    cd(path);
    copyfile(strcat(p,"T1_nyul.nii"),fullfile(myelin_path,char(p),'MNI_T1nmlz.nii'))
    path = normalizedpath2;
    cd(path);
    copyfile(strcat(p,"T2_nyul.nii"),fullfile(myelin_path,char(p),'MNI_T2nmlz.nii'))
end
%%
subjID=["XXXX"];
normalizedpath1 = 'Z:\Imaging\Multimodal\Myelin\nyul_normalizedT1';
normalizedpath2 = 'Z:\Imaging\Multimodal\Myelin\nyul_normalizedT2';
myelin_path = 'Z:\Imaging\Multimodal\Myelin\Normal';
for p = subjID
    path = normalizedpath1;
    cd(path);
    copyfile(strcat(p,"T1_nyul.nii"),fullfile(myelin_path,char(p),'MNI_T1nmlz.nii'))
    path = normalizedpath2;
    cd(path);
    copyfile(strcat(p,"T2_nyul.nii"),fullfile(myelin_path,char(p),'MNI_T2nmlz.nii'))
end
%% regenerate T1oT2files
subjID=["XXXX"];
myelin_path = 'Z:\Imaging\Multimodal\Myelin\Normal';
load('dummy_header_MRF_1_0mm.mat');
for p = subjID
    path = strcat(myelin_path,'\',p);
    cd(path);
    T1_raw = load_untouch_nii('MNI_T1nmlz.nii');
    T2_raw = load_untouch_nii('MNI_T2nmlz.nii');

    img1 = single(T1_raw.img);
    mask = single((img1 > mean(img1(:))));
    img3 = single(T2_raw.img);
    t1overt2 = (img1./img3).*mask;
    t1overt2(isnan(t1overt2(:)))= 0;
    t1overt2(isinf(t1overt2(:)))= 0;
   
    finalimage = dummy_1_0mm_un;
    finalimage.img = t1overt2;
%     finalimage.fileprefix ='T1overT2';
%     finalimage.hdr.dime.datatype = 16;
%     finalimage.hdr.dime.bitpix = 32;
    save_untouch_nii(finalimage, 'MNI_T1oT2nmlz.nii');
end
%%
subjID=["XXX"];
myelin_path = 'Z:\Imaging\Multimodal\Myelin\Patients';
for p = subjID
    path = strcat(myelin_path,'\',p);
    cd(path);
    copyfile MNI_T1im.nii MNI_T1oT2nmlz.nii
    T1_raw = load_untouch_nii('MNI_T1nmlz.nii');
    T2_raw = load_untouch_nii('MNI_T2nmlz.nii');
    finalimage= load_untouch_nii('MNI_T1oT2nmlz.nii');

    img1 = double(T1_raw.img);
    mask = double((img1 > mean(img1(:))));
    img3 = double(T2_raw.img);
    t1overt2 = (img1./img3).*mask;
    t1overt2(isnan(t1overt2(:)))= 0;
    t1overt2(isinf(t1overt2(:)))= 0;
    finalimage.img = t1overt2;
    finalimage.fileprefix ='T1overT2';
    finalimage.hdr.dime.datatype = 16;
    finalimage.hdr.dime.bitpix = 32;
    save_untouch_nii(finalimage, 'MNI_T1oT2nmlz.nii');
end
%%

folderPath = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\postop_analysis';
items = dir(folderPath);
folders = items([items.isdir]);
folderNames = {folders.name};
folderNames = folderNames(~ismember(folderNames, {'.', '..'}));
folderMatrix = string(folderNames); 
folderPath = 'C:\postopanalysis'
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
for p = folderMatrix
    cd(folderPath)
    mkdir(p)
    path = strcat(folderPath,'\',p);
    cd(path);
    mkdir Pre
    mkdir Post
    copyfile(strcat(MRF_path,'\',p,'\VBM_newANN\T1_cat.nii'), strcat(path,'\Pre\MPRAGE.nii'));
    copyfile(strcat(MRF_path,'\',p,'\postop.nii'), strcat(path,'\Post\postop.nii'));
end
%% dilate images
folderPath = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\postop_analysis';
items = dir(folderPath);
folders = items([items.isdir]);
folderNames = {folders.name};
folderNames = folderNames(~ismember(folderNames, {'.', '..'}));
folderMatrix = string(folderNames); 
folderPath = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients'
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
for p = folderMatrix
    cd(strcat(MRF_path,'\',p,'\MRF_VBM'))
    copyfile MNI_resected_binary.nii MNI_resected_dilated.nii
    file = load_untouch_nii('MNI_resected_dilated.nii');
    image = single(file.img);
    se = strel('sphere',5);
    dilated = imdilate(image,se);
    file.img = dilated;
    save_untouch_nii(file, 'MNI_resected_dilated.nii')
end
%% 
subjID=["XXX"];
subjID=["XXX"];
MRF_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
train_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\imagesTr';
label_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\labelsTr';
unziptrain_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\imagesTr_unzipped';
unziplabel_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\labelsTr_unzipped';
for p = subjID
    cd(strcat(MRF_path,'\',p,'\MRF_VBM'));
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_T1w.nii'), strcat(unziptrain_path,'\', p,'_0000.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0000.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_T1.nii'), strcat(unziptrain_path,'\', p,'_0001.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0001.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_T2.nii'), strcat(unziptrain_path,'\', p,'_0002.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0002.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_GM_prob.nii'), strcat(unziptrain_path,'\', p,'_0003.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0003.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_WM_prob.nii'), strcat(unziptrain_path,'\', p,'_0004.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0004.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_CSF_prob.nii'), strcat(unziptrain_path,'\', p,'_0005.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0005.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_ROI_final.nii'), strcat(unziplabel_path,'\', p,'.nii'))
    gzip(strcat(unziplabel_path,'\', p,'.nii'), label_path)
end


%%
subjID=["XXXX"];
% subjID=["P01_12251"];
MRF_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
train_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\imagesTr';
label_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\labelsTr';
unziptrain_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\imagesTr_unzipped';
unziplabel_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\labelsTr_unzipped';
for p = subjID
    cd(strcat(MRF_path,'\',p,'\MRF_VBM'));
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_T1w.nii'), strcat(unziptrain_path,'\', p,'_0000.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0000.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_T1.nii'), strcat(unziptrain_path,'\', p,'_0001.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0001.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_T2.nii'), strcat(unziptrain_path,'\', p,'_0002.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0002.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_GM_prob.nii'), strcat(unziptrain_path,'\', p,'_0003.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0003.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_WM_prob.nii'), strcat(unziptrain_path,'\', p,'_0004.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0004.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_CSF_prob.nii'), strcat(unziptrain_path,'\', p,'_0005.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0005.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_ROI_final.nii'), strcat(unziplabel_path,'\', p,'.nii'))
    gzip(strcat(unziplabel_path,'\', p,'.nii'), label_path)
end


%%
subjID=["XXXX"];
MRF_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
% train_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset013_MRFzclin\imagesTr';
% label_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset013_MRFzclin\labelsTr';
train_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset002_MRF_small\imagesTr';
label_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset002_MRF_small\labelsTr';
unziptrain_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\imagesTr_unzipped';
unziplabel_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\labelsTr_unzipped';
for p = subjID
    cd(strcat(MRF_path,'\',p,'\MRF_VBM'));
    p
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_T1w.nii'), strcat(unziptrain_path,'\', p,'_0000.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0000.nii'), train_path)    
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_T1.nii'), strcat(unziptrain_path,'\', p,'_0001.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0001.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_T2.nii'), strcat(unziptrain_path,'\', p,'_0002.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0002.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_GM_prob.nii'), strcat(unziptrain_path,'\', p,'_0003.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0003.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_WM_prob.nii'), strcat(unziptrain_path,'\', p,'_0004.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0004.nii'), train_path)    
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_CSF_prob.nii'), strcat(unziptrain_path,'\', p,'_0005.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0005.nii'), train_path)
%     copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_WM_prob.nii'), strcat(unziptrain_path,'\', p,'_0006.nii'))
%     gzip(strcat(unziptrain_path,'\', p,'_0006.nii'), train_path)
%     copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_CSF_prob.nii'), strcat(unziptrain_path,'\', p,'_0007.nii'))
%     gzip(strcat(unziptrain_path,'\', p,'_0007.nii'), train_path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_ROI_final.nii'), strcat(unziplabel_path,'\', p,'.nii'))
    gzip(strcat(unziplabel_path,'\', p,'.nii'), label_path)

end

%% for 6 fold validations copy files to imagetest set
subjID=["XXX"];
MRF_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
train_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset057_6fold32\imagesTr';
test_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset057_6fold32\imagesTs';
label_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset057_6fold32\labelsTr';
label_test = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset057_6fold32\labelsTs';

for p = subjID
    cd(strcat(MRF_path,'\',p,'\MRF_VBM'));
    p
    movefile(strcat(train_path,'\',p,'_0000.nii.gz'), test_path)
    movefile(strcat(train_path,'\',p,'_0001.nii.gz'), test_path)
    movefile(strcat(train_path,'\',p,'_0002.nii.gz'), test_path)
    movefile(strcat(train_path,'\',p,'_0003.nii.gz'), test_path)
    movefile(strcat(train_path,'\',p,'_0004.nii.gz'), test_path)
    movefile(strcat(train_path,'\',p,'_0005.nii.gz'), test_path)
    movefile(strcat(train_path,'\',p,'_0006.nii.gz'), test_path)
    movefile(strcat(train_path,'\',p,'_0007.nii.gz'), test_path)
    movefile(strcat(label_path,'\',p,'.nii.gz'), label_test)
end

%%
% subjID=["XXXX"];
% subjID=["P01_12251"];
MRF_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
train_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset009_MRF_clin\imagesTr';
label_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\labelsTr';
unziptrain_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\imagesTr_unzipped';
unziplabel_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\labelsTr_unzipped';
for p = subjID
    cd(strcat(MRF_path,'\',p,'\MRF_VBM'));
    p
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_MPRAGE_fn.nii'), strcat(unziptrain_path,'\', p,'_0000.nii'))
    gzip(strcat(unziptrain_path,'\', p,'_0000.nii'), train_path)   
%     copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_T1.nii'), strcat(unziptrain_path,'\', p,'_0001.nii'))
%     gzip(strcat(unziptrain_path,'\', p,'_0001.nii'), train_path)
%     copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_T2.nii'), strcat(unziptrain_path,'\', p,'_0002.nii'))
%     gzip(strcat(unziptrain_path,'\', p,'_0002.nii'), train_path)
%     copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_age.nii'), strcat(unziptrain_path,'\', p,'_0003.nii'))
%     gzip(strcat(unziptrain_path,'\', p,'_0003.nii'), train_path)
%     copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_gender.nii'), strcat(unziptrain_path,'\', p,'_0004.nii'))
%     gzip(strcat(unziptrain_path,'\', p,'_0004.nii'), train_path)    
%     copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_GM_prob.nii'), strcat(unziptrain_path,'\', p,'_0005.nii'))
%     gzip(strcat(unziptrain_path,'\', p,'_0005.nii'), train_path)
%     copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_WM_prob.nii'), strcat(unziptrain_path,'\', p,'_0006.nii'))
%     gzip(strcat(unziptrain_path,'\', p,'_0006.nii'), train_path)
%     copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_CSF_prob.nii'), strcat(unziptrain_path,'\', p,'_0007.nii'))
%     gzip(strcat(unziptrain_path,'\', p,'_0007.nii'), train_path)
end
%%
subjID=["XXX"];
MRF_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\';
train_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset015_MRFz\imagesTr';
label_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset015_MRFz\labelsTr';
unziptrain_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\imagesTr_unzipped';
unziplabel_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\labelsTr_unzipped';
for p = subjID
    cd(strcat(MRF_path,'\',p,'\MRF_VBM'));
    p
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_MPRAGE_fn.nii'), strcat(MRF_path, p, '\MRF_VBM\MNI_MPRAGE_brain.nii'))
    T1 = load_untouch_nii('MNI_MPRAGE_brain.nii');
    T1im = double(T1.img);
    mask = double(load_untouch_nii('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\MNI152_T1_1mm_brain.nii').img);
    mask(mask>0) = 1;
    T1im = T1im .* mask;
    T1.img = T1im;
    save_untouch_nii(T1, 'MNI_MPRAGE_brain.nii');

end
%% check nan
subjID=["XXX"];
MRF_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\';
train_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset015_MRFz\imagesTr';
label_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset015_MRFz\labelsTr';
unziptrain_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\imagesTr_unzipped';
unziplabel_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\labelsTr_unzipped';
for p = subjID
    cd(strcat(MRF_path,'\',p,'\MRF_VBM'));
    p
    T1 = load_untouch_nii('MNI_T1LMz.nii');
    T1im = double(T1.img);
    sum(isnan(T1im), 'all')

end
%%
subjID=["XXXX"];
MRF_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
% train_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset013_MRFzclin\imagesTr';
% label_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset013_MRFzclin\labelsTr';
train_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset016_MRFz\imagesTr';
label_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\nnUNet_raw\Dataset016_MRFz\labelsTr';
unziptrain_path = 'Z:\Imaging\Multimodal\MRF\Peter\nnUNet_CV\Dataset016_MRFz\valNN_probs';
unziplabel_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\MRF_to_nnunet\labelsTr_unzipped';
for p = subjID
    cd(strcat(MRF_path,'\',p,'\MRF_VBM'));
    p
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_lesionprob_ANN_fn.nii'), strcat(unziptrain_path,'\', p,'_valnn.nii'))
end
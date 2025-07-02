subjID=[XXXX];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';

for p = subjID
    p
    path = strcat(MRF_path,'\',p);
    cd(path)
    mkdir(fullfile(path,'VBM_newANN'));
    dest = fullfile(path,'VBM_newANN');
    cd('VBM')
    sourceFile = 'T1.nii';
    destFile = fullfile(dest, 'T1_cat.nii'); 
    copyfile(sourceFile, destFile);
end





%% perform VBM with new ANN
subjID=[XXXX];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';

for p = subjID
    p
    path = strcat(MRF_path,'\',p);
    dest = fullfile(path,'VBM_newANN');
    cd(dest)
    map18('mapFCD','chunkedNN');
end


%%
subjID=[XXXX];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

for p = subjID
    p
    path = strcat(MRF_path,'\',p);
    cd(path)
    mkdir(fullfile(path,'VBM_newANN'));
    dest = fullfile(path,'VBM_newANN');
    cd('reorient_all')
    sourceFile = 'T1.nii';
    destFile = fullfile(dest, 'T1_cat.nii'); 
    copyfile(sourceFile, destFile);
    
end

%% perform VBM with new ANN

subjID=[XXX];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

for p = subjID
    path = strcat(MRF_path,'\',p);
    dest = fullfile(path,'VBM_newANN');
    cd(dest)
    map18('mapFCD','valNN');
    map18('inverse');
end

%%
subjID=[XXX];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
for p = subjID
    p
    path = strcat(MRF_path,'\',p);
    dest = fullfile(path,'VBM_newANN');
    cd(dest)
    map18('mapFCD','valNN');
    map18('inverse');
end
%% warp back to native space
subjID=[XXX];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
for p = subjID
    path = strcat(MRF_path,'\',p);
    dest = fullfile(path,'VBM_newANN');
    cd(dest)
    map18('inverse');
end

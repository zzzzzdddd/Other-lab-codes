%% Code for IMA to Nifti with adding header information (12/14/2021; J)
clear all
%% T1 NEW MRF
subjID = [XXXX];
scode = [XXX];
images = ["_FLAIRmaps.IMA" "_GMmaps.IMA" "_T1maps.IMA" "_T1W.IMA" "_T2maps.IMA" "_WMFmaps.IMA" "_WMmaps.IMA"];
nf = ["FLAIR_NEW_data.nii" "GM_NEW_data.nii" "T1_NEW_data.nii" "T1W_NEW_data.nii" "T2_NEW_data.nii" "WMF_NEW_data.nii" "WM_NEW_data.nii"];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';

for v = subjID
    path = strcat(MRF_path,'\',v,'\fastMRF');
    cd(path);
    for i = images
        f = char(strcat(scode(find(subjID==v)),i));
        T1 = dicomread(f);
        T1 = single(permute(T1,[1 2 4 3]));
        T1 = squeeze(T1);
        T1 = flip(T1,3);
        T1 = rot90(T1,-1);
        
        newfilename = char(nf(find(images==i)));
        fid = fopen(newfilename,'w');
        fwrite(fid, T1, 'float');
        fclose(fid);
        nii = make_nii(T1);

        load('dummy_header_MRF_1_0mm.mat');
        
        final_header = dummy_1_0mm_un;
        final_header.img = nii.img;
        save_untouch_nii(final_header, newfilename);
    end
end

%% T1 OLD MRF
subjID = [XXX];
scode = [XXXX];
images = ["_480_old_T1maps.IMA" "_480_old_T1W.IMA" "_480_old_T2maps.IMA" "_1200_old_T1maps.IMA" "_1200_old_T1W.IMA" "_1200_old_T2maps.IMA"];
nf = ["T1_OLD_480_data.nii" "T1W_OLD_480_data.nii" "T2_OLD_480_data.nii" "T1_OLD_1200_data.nii" "T1W_OLD_1200_data.nii" "T2_OLD_1200_data.nii"];
for v = subjID
    path = strcat(MRF_path,'\',v,'\fastMRF');
    cd(path);
    for i = images
        f = char(strcat(scode(find(subjID==v)),i));
        if exist(f)
            T1 = dicomread(f);
            T1 = single(permute(T1,[1 2 4 3]));
            T1 = squeeze(T1);
            T1 = flip(T1,3);
            T1 = rot90(T1,-1);
            
            newfilename = char(nf(find(images==i)));
            fid = fopen(newfilename,'w');
            fwrite(fid, T1, 'float');
            fclose(fid);
            nii = make_nii(T1);
            
            load('dummy_header_MRF_1_0mm.mat');
            
            final_header = dummy_1_0mm_un;
            final_header.img = nii.img;
            save_untouch_nii(final_header, newfilename);
        end
    end
end


subjID=[XXXX];

MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
cd('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T');
for p = subjID
    p
    path = strcat(MRF_path,'\',p);
    cd(path)
    copyfile T1w_data.nii motion_mask.nii
    binaryimage = load_untouch_nii('motion_mask.nii');
    final_area = double(binaryimage.img);
    final_area(:,:,:) = 1;
    switch p
        case "XXXX"
            final_area(:,:,:) = 0;
        case "XXXX"
            final_area(:,:,:) = 0;
        case "XXXX"
            final_area(:,:,:) = 0;
        case "XXXX"
            final_area(:,:,:) = 0;
        case "XXXX"
            final_area(:,:,38:77) = 0;        
        case "XXXX"
            final_area(:,:,:) = 0;
        case "XXX"
            final_area(:,:,85:end) = 0;
        case "XXXX"
            final_area(:,:,:) = 0;
        case "XXXX"
            final_area(:,:,1:52) = 0;            
        case "XXXX"
            final_area(:,:,128:129) = 0;
        case "XXXXX"
            final_area(:,:,114:115) = 0;
            final_area(:,:,43:50) = 0;
        case "XXXX"
            final_area(:,:,41:45) = 0;
            final_area(:,:,111:112) = 0;
        case "XXXX"
            final_area(:,:,101:102) = 0;
        case "XXXXX"
            final_area(:,:,79:80) = 0;
        case "XXX"
            final_area(:,:,64:77) = 0;
    end
    binaryimage.img = final_area;
    save_untouch_nii(binaryimage, 'motion_mask.nii')
end
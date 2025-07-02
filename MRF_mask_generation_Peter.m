%% T1 + T1w  (old)   480

mask1brain = load_untouch_nii('T1_data_brain_ini_mask.nii'); 'T1_OLD_480_data.nii'
mask1_img = single(mask1brain.img);

mask2brain = load_untouch_nii('T1w_data_brain_ini_mask.nii'); 
mask2_img = single(mask2brain.img);

T1_mask_er = erode_3d(mask1_img, 4);
[x y z] = size(T1_mask_er);

combine_mask = zeros(x,y,z); 
combine_mask(:,y/2:end,:) = mask2_img(:,y/2:end,:);
combine_mask(:,1:y/2-1,:) = T1_mask_er(:,1:y/2-1,:);

fid = fopen('MRF_brain_mask.nii','w'); %% This
fwrite(fid, combine_mask, 'float'); 
fclose(fid);
mask1brain.img = combine_mask;
save_untouch_nii(mask1brain, 'MRF_brain_mask.nii');


%% For new MRF T1w two different f values
subjID = [XXXX];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';

for v = subjID
    path = strcat(MRF_path,'\',v,'\fastMRF');
    cd(path);
    
    mask1brain = load_untouch_nii('T1w_NEW_data_0.25_brain_ini.nii');  
    brain1_img = single(mask1brain.img);
    brain1_img(brain1_img>=80) = 0;

    mask2brain = load_untouch_nii('T1w_NEW_data_0.3_brain_ini.nii');  
    brain2_img = single(mask2brain.img);
    brain2_img(brain2_img>=80) = 0;

    mask1 = load_untouch_nii('T1w_NEW_data_0.25_brain_ini_mask.nii');  
    mask1_img = single(mask1.img);
    mask1_img = mask1_img.*brain1_img;
    mask1_img(mask1_img>0) = 1;

    mask2 = load_untouch_nii('T1w_NEW_data_0.3_brain_ini_mask.nii');  
    mask2_img = single(mask2.img);
    mask2_img = mask2_img.*brain1_img;
    mask2_img(mask2_img>0) = 1;
    
%     se = strel('sphere',2);
%     mask1_img = imerode(mask1_img, se);
%     mask2_img = imerode(mask2_img, se);
% 
%     mask1_img = ExtractNLargestBlobs(mask1_img, 1);
%     mask2_img = ExtractNLargestBlobs(mask2_img, 1);
% 
%     mask1_img = imdilate(mask1_img, se);
%     mask2_img = imdilate(mask2_img, se);

    [x y z] = size(mask1_img);
    
    combine_mask = zeros(x,y,z); 
    combine_mask(:,:,1:50) = mask2_img(:,:,1:50);
    combine_mask(:,:,51:end) = mask1_img(:,:,51:end);

    se = strel("disk",3)
    for i = 1:z
        combine_mask(:,:,i) = imerode(combine_mask(:,:,i),se);
        combine_mask(:,:,i) = ExtractNLargestBlobs2d(combine_mask(:,:,i), 1);
        combine_mask(:,:,i) = imdilate(combine_mask(:,:,i),se);
    end
 
    fid = fopen('MRF_brain_mask_new.nii','w'); 
    fwrite(fid, combine_mask, 'float');
    fclose(fid);
    mask1brain.img = combine_mask;
    save_untouch_nii(mask1brain, 'MRF_brain_mask_new.nii');
end

function binaryImage = ExtractNLargestBlobs(binaryImage, numberToExtract)
	[labeledImage, numberOfBlobs] = bwlabeln(binaryImage);
	blobMeasurements = regionprops3(labeledImage, 'Volume');
	% Get all the areas
	allAreas = [blobMeasurements.Volume];
	if numberToExtract > 0
		[sortedAreas, sortIndexes] = sort(allAreas, 'descend');
	end
	% Extract the "numberToExtract" largest blob(a)s using ismember().
	biggestBlob = ismember(labeledImage, sortIndexes(1:numberToExtract));
	% Convert from integer labeled image into binary (logical) image.
	binaryImage = biggestBlob > 0;
end

function binaryImage = ExtractNLargestBlobs2d(binaryImage, numberToExtract)
	[labeledImage, numberOfBlobs] = bwlabeln(binaryImage);
	blobMeasurements = regionprops(labeledImage, 'Area');
	% Get all the areas
	allAreas = [blobMeasurements.Area];
	if numberToExtract > 0
		[sortedAreas, sortIndexes] = sort(allAreas, 'descend');
    end
    s = size(sortedAreas);
	if s(1) == 0
        biggestBlob = binaryImage;
    else
        biggestBlob = ismember(labeledImage, sortIndexes(1:numberToExtract));
        binaryImage = biggestBlob > 0;
    end
%     biggestBlob = ismember(labeledImage, sortIndexes(1:numberToExtract));
%     binaryImage = biggestBlob > 0;
end

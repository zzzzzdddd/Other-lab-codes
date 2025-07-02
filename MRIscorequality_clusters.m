% clear all
subjID=[XXX
     ];

MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
% MRF_path= '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Patients';
% subjID=[XXXX];
% subjID=[XXXX];
% 
% MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';

for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni = single(ann.img);
    anni(anni>0) = 1;
    mask = single(load_untouch_nii('MNI_Brain_Mask.nii').img);
    mask(mask>0.95) = 1;
    copyfile MNI_ROI_final.nii MNI_cvGM_clusters.nii
    copyfile MNI_ROI_final.nii MNI_cvWM_clusters.nii
    copyfile MNI_ROI_final.nii MNI_cvcGMWM_clusters.nii

    anni = anni.*mask;
    GM = single(load_untouch_nii('MNI_GM_prob.nii').img);
    GM(GM<0.8) = 0;
    GM(GM>=0.8) = 1;
    WM = single(load_untouch_nii('MNI_WM_prob.nii').img);
    WM(WM<0.8) = 0;
    WM(WM>=0.8) = 1;

    MRI = load_untouch_nii('MNI_MPRAGE.nii');
    image = single(MRI.img);

    intensitiesGM = image.*GM;
    intensitiesWM = image.*WM;

    [labeledImage, numberOfBlobs] = bwlabeln(anni.*GM);
	blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
	allinds = [blobMeasurements.VoxelIdxList];

%     roi = load_untouch_nii('MNI_ROI_final.nii');
%     roii = single(roi.img);
%     [roiblob, nBlobs] = bwlabeln(roii);
%     roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
%     roiinds = roiMeasurements.VoxelIdxList{1};
    
    trueind = [];
%     for i = 1:numberOfBlobs
%         A = blobMeasurements.VoxelIdxList{i};
%         if sum(ismember(roiinds,A))>0
%             trueind = horzcat(trueind, i);
%         end
%     end
    for i = 1:numberOfBlobs
%         if (size(allinds{i}) >= 10)
%         if (ismember(i,trueind) == 0) & (size(allinds{i}) >= 10)
            meanGM = mean(intensitiesGM(allinds{i}));
            stdGM = std(intensitiesGM(allinds{i}));
            cvGM  = stdGM./meanGM; % coefficient of variation for GM signal intensities

            file = load_untouch_nii('MNI_cvGM_clusters.nii');
            file.img(allinds{i}) = cvGM;
            save_untouch_nii(file, 'MNI_cvGM_clusters.nii');
%         end
    end

    [labeledImage, numberOfBlobs] = bwlabeln(anni.*WM);
	blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
	allinds = [blobMeasurements.VoxelIdxList];

%     roi = load_untouch_nii('MNI_ROI_final.nii');
%     roii = single(roi.img);
%     [roiblob, nBlobs] = bwlabeln(roii);
%     roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
%     roiinds = roiMeasurements.VoxelIdxList{1};
    
    trueind = [];
%     for i = 1:numberOfBlobs
%         A = blobMeasurements.VoxelIdxList{i};
%         if sum(ismember(roiinds,A))>0
%             trueind = horzcat(trueind, i);
%         end
%     end
    for i = 1:numberOfBlobs
%         if (size(allinds{i}) >= 10)
%         if (ismember(i,trueind) == 0) & (size(allinds{i}) >= 10)
            meanWM = mean(intensitiesWM(allinds{i}));
            stdWM = std(intensitiesWM(allinds{i}));
            cvWM       = stdWM./meanWM; % coefficient of variation for WM signal intensities

            file = load_untouch_nii('MNI_cvWM_clusters.nii');
            file.img(allinds{i}) = cvWM;
            save_untouch_nii(file, 'MNI_cvWM_clusters.nii');


%         end
    end

    [labeledImage, numberOfBlobs] = bwlabeln(anni.*mask);
	blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
	allinds = [blobMeasurements.VoxelIdxList];

%     roi = load_untouch_nii('MNI_ROI_final.nii');
%     roii = single(roi.img);
%     [roiblob, nBlobs] = bwlabeln(roii);
%     roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
%     roiinds = roiMeasurements.VoxelIdxList{1};
    
    trueind = [];
%     for i = 1:numberOfBlobs
%         A = blobMeasurements.VoxelIdxList{i};
%         if sum(ismember(roiinds,A))>0
%             trueind = horzcat(trueind, i);
%         end
%     end
    for i = 1:numberOfBlobs
%         if (size(allinds{i}) >= 10)
%         if (ismember(i,trueind) == 0) & (size(allinds{i}) >= 10)
            meanWM = mean(intensitiesWM(allinds{i}));
            stdWM = std(intensitiesWM(allinds{i}));
            cvWM       = stdWM./meanWM; % coefficient of variation for WM signal intensities
            meanGM = mean(intensitiesGM(allinds{i}));
            stdGM = std(intensitiesGM(allinds{i}));
            cvGM       = stdGM./meanGM; % coefficient of variation for GM signal intensities
            cvcGMWM    = abs(meanGM-meanWM)./mean([stdGM stdWM]);
            file = load_untouch_nii('MNI_cvcGMWM_clusters.nii');
            file.img(allinds{i}) = cvcGMWM;
            save_untouch_nii(file, 'MNI_cvcGMWM_clusters.nii');
%         end
    end
        
end
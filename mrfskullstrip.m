subjID=[XXXX];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

% subjID=[XXXXX];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
for p = subjID
%     path = strcat(MRF_path,'\',p);
%     path = strcat(MRF_path,'\',p,'\MRF_new');
    path = strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\TJU-Jefferson\',p,'\MRF_new');
    cd(path)
    gunzip('m0.nii.gz')
    m0 = double(load_untouch_nii('m0.nii').img);
    m0(m0>0)=1;
    binaryimage = load_untouch_nii('m0_mask.nii');
    mpg = single(load_untouch_nii('T1w_data.nii').img);
    maskrevise = m0;
    maskrevise(mpg>=0.8)=0;
    maskrevise(mpg==0)=0;
    for z = 1:144
        maskrevise(:,:,z) = bwareafilt(logical(maskrevise(:,:,z)), 3, 4);        
    end
    maskrevise = imfill(maskrevise, 'holes');
    binaryimage.img = maskrevise;
    save_untouch_nii(binaryimage,'m0_mask.nii');
end

%%
copyfile T1w_data.nii skullbound.nii
thresh = 0.01
BW = edge3(t1wi,"approxcanny",thresh);
skb.img = BW;
save_untouch_nii(skb,'skullbound.nii');
%%
copyfile T1w_data.nii skullmask.nii
t1w = load_untouch_nii('T1w_data.nii');
t1wi = double(t1w.img);
% t1s = load_untouch_nii('T1_data.nii');
% t1si = double(t1s.img);
% t2s = load_untouch_nii('T2_data.nii');
% t2si = double(t2s.img);
sk = load_untouch_nii('skullmask.nii');
ski = double(sk.img);
ski = ski - BW;
ski(ski<0.4) = 0;
ski((wmi+gmi+csfi)>0.95) = 0;
[labeledImage, numberOfBlobs] = bwlabeln(ski);
blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
% allinds = [blobMeasurements.VoxelIdxList(blobMeasurements.Volume>10000)];
allinds = [blobMeasurements.VoxelIdxList(blobMeasurements.Volume==max(blobMeasurements.Volume))];
ski(allinds{1}) = 1;
ski(ski<1)=0;
ski = ski + BW;
ski(ski>1)=1;

se = strel('sphere',1);
imerode(ski,se);
imerode(ski,se);
imdilate(ski,se);
imdilate(ski,se);
sk.img = ski;
save_untouch_nii(sk,'skullmask.nii');
%%
[labeledImage, numberOfBlobs] = bwlabeln(ski);
blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
allinds = [blobMeasurements.VoxelIdxList(blobMeasurements.Volume==max(blobMeasurements.Volume))];
ski(allinds{1}) = 1;
ski(ski<1)=0;
sk.img = ski;
save_untouch_nii(sk,'skullmask.nii');
%% start over
copyfile T1w_data.nii skullmask.nii

t1max = 1700; 
t2max = 180;
t1min = 750;
t2min = 30;
sk = load_untouch_nii('skullmask.nii');
ski = double(sk.img);
ski((t1si<t1max)&(t1si>t1min)&(t2si<t2max)&(t2si>t2min)) = 0;
ski((ski<0.3)&(ski>0.2)) = 0;
% ski((wmi+gmi+csfi)>0.95) = 0;
ski(ski>0) = 1;

ski = ski - BW;
[labeledImage, numberOfBlobs] = bwlabeln(ski);
blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
% allinds = [blobMeasurements.VoxelIdxList(blobMeasurements.Volume>10000)];
allinds = [blobMeasurements.VoxelIdxList(blobMeasurements.Volume==max(blobMeasurements.Volume))];
ski(allinds{1}) = 1;
ski = ski + BW;
ski(ski>1)=1;
sk.img = ski;
save_untouch_nii(sk,'skullmask.nii');

%% WM keep middle
c1 = load_untouch_nii('c1T1w_data.nii');
c1i = double(c1.img);
c2 = load_untouch_nii('c2T1w_data.nii');
c2i = double(c2.img);
c3 = load_untouch_nii('c3T1w_data.nii');
c3i = double(c3.img);
copyfile T1w_data.nii skullmask.nii
sk = load_untouch_nii('skullmask.nii');
ski = double(sk.img);
ski = c1i+c2i+c3i;
ski(ski>0.98)=1;
ski(ski<1)=0;
sk.img = ski;
save_untouch_nii(sk,'skullmask.nii');

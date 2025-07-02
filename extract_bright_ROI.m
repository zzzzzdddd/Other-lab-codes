clear all
data = load_untouch_nii('overlay_MNI_improve.nii');
im = single(data.img);

fid = fopen('ROI_MNI_improve1.nii','w'); 

[mxv,idx] = max(im(:));
level = graythresh(im)*0.8;
bina = zeros(size(im));
bina(find(im>(level*mxv))) = 1;
bina = single(bwmorph3(bina,'clean'));
niftiwrite(bina,'ROI_MNI_improve1.nii')

fid2 = fopen('ROI_MNI_improve2.nii','w'); 
SE = strel('disk',2);
bina = imclose(bina,SE);
bina = imopen(bina,SE);
niftiwrite(bina,'ROI_MNI_improve2.nii')
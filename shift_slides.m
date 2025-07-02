clear all
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_new\
fileList = dir('*.nii');
files = string({fileList.name});
topslice = 46;

for f = files
    structure = load_untouch_nii(char(f));
    copied = single(structure.img);
    totslice = size(copied,3);
    structure.img(:,:,1:(totslice-topslice+1)) = copied(:,:,topslice:end);
    structure.img(:,:,(totslice-topslice+1):end) = copied(:,:,1:topslice);
    save_untouch_nii(structure,char(f));
end
%%
clear all
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXXX\MRF_new
fileList = dir('*.nii');
files = string({fileList.name});
topslice = 44;

for f = files
    structure = load_untouch_nii(char(f));
    copied = single(structure.img);
    totslice = size(copied,3);
    structure.img(:,:,1:(totslice-topslice+1)) = copied(:,:,topslice:end);
    structure.img(:,:,(totslice-topslice+1):end) = copied(:,:,1:topslice);
    save_untouch_nii(structure,char(f));
end
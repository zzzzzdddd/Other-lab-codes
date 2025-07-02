subjID=[XXXX];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';

for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni =  single(ann.img);
    anni(anni>0) = 1;
    anni = -(anni - 1);
    cd('C:\Users\Irene\Downloads\MRF_VBM\imagesTr');
    file = fullfile('C:\Users\Irene\Downloads\MRF_VBM\imagesTr',strcat(p,".nii"))
    im = load_untouch_nii(char(file));
    proim = single(im.img);
    for l = 1:(size(proim,4))
        proim(:,:,:,l) = proim(:,:,:,l).*anni;
    end 
    im.img = proim;
    save_untouch_nii(im, file)
end

%%
cd C:\Users\Irene\Downloads
FileData = load('2Dtraining.mat');
csvwrite('traindata.csv', FileData.xtrain,FileData.re);
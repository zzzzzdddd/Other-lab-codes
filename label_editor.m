
labels = [XXX];
ROIs = [XXXX]
for i = 1:size(labels,2)
    tf = char(fullfile('C:\Users\Irene\Downloads\MRF_VBM\imagesTs',labels(i)));
    rf = char(fullfile('C:\Users\Irene\Downloads\MRF_VBM\imagesTs',ROIs(i)));
    dat = load_untouch_nii(tf);
    dat2 = load_untouch_nii(rf);
    dat.img = dat.img .* dat2.img;
    save_untouch_nii(dat,tf)
end
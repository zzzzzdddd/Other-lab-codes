clear all

%% MNI space ROI
cd T:\Imaging\VBM\Phelps_Molly_89152930_20211225_20031120

junc = load_untouch_nii('wT1_junction_z_score.nii');
juncz = single(junc.img);

ext = load_untouch_nii('wT1_extension_z_score.nii');
extz = single(ext.img);

thick = load_untouch_nii('wT1_thickness_z_score.nii');
thickz = single(thick.img);

juncROI = zeros(size(juncz));
extROI= zeros(size(juncz));
thickROI = zeros(size(juncz));
combineROI = zeros(size(juncz));
combineROI_strict = zeros(size(juncz));

juncROI(juncz>=4) = 1;
extROI(extz>=6) = 1;
thickROI(thickz>=4) = 1;
combineROI((juncROI==1)|(extROI==1)|(thickROI==1)) = 1;
combineROI_strict(((juncROI==1)&(extROI==1))|((extROI==1)&(thickROI==1))|((juncROI==1)&(thickROI==1))) = 1;

niftiwrite(juncROI, 'MNIjunctionROI.nii');
niftiwrite(extROI, 'MNIextROI.nii');
niftiwrite(thickROI, 'MNIthickROI.nii');
niftiwrite(combineROI,'MNIcombineROI.nii');
niftiwrite(combineROI_strict,'MNIcombineROI_strict.nii');

%% original space
junc = load_untouch_nii('T1_junction_z_score.nii');
juncz = single(junc.img);

ext = load_untouch_nii('T1_extension_z_score.nii');
extz = single(ext.img);

thick = load_untouch_nii('T1_thickness_z_score.nii');
thickz = single(thick.img);

juncROI = zeros(size(juncz));
extROI= zeros(size(juncz));
thickROI = zeros(size(juncz));
combineROI = zeros(size(juncz));
combineROI_strict = zeros(size(juncz));

juncROI(juncz>=4) = 1;
extROI(extz>=6) = 1;
thickROI(thickz>=4) = 1;
combineROI((juncROI==1)|(extROI==1)|(thickROI==1)) = 1;
combineROI_strict(((juncROI==1)&(extROI==1))|((extROI==1)&(thickROI==1))|((juncROI==1)&(thickROI==1))) = 1;

niftiwrite(juncROI, 'junctionROI.nii');
niftiwrite(extROI, 'extROI.nii');
niftiwrite(thickROI, 'thickROI.nii');
niftiwrite(combineROI,'combineROI.nii');
niftiwrite(combineROI_strict,'combineROI_strict.nii');

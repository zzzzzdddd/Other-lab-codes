clear all
% 
subjID=[XXXX];

MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
for p = subjID
    path1 = strcat(MRF_path,'\',p,'\MRF_VBM');
    path = strcat(MRF_path,'\',p,'\MRFz');
    if exist(path)
       cd(path);
       copyfile('MNI_T1.nii', path1)
       copyfile('MNI_T2.nii', path1)
       copyfile('MNI_GM_prob.nii', path1)
       copyfile('MNI_WM_prob.nii', path1)
    end
end

%% unzip further
subjID=[XXX];

MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
for p = subjID
    path = strcat(MRF_path,'\',p,'\VBM');
    if exist(path)
       cd(path);
       unzip('T1_further_morphometric_results.zip');
    end
end

%% store all T1 data and get mean and std

subjID=[XXX];

MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
mT1 = zeros(182,218,182);
sdT1 = mT1;
mT2 = mT1;
sdT2 = mT1;
mGM = mT1;
sdGM = mT1;
mWM = mT1;
sdWM = mT1;
n = 1;
for p = subjID
    path = strcat(MRF_path,'\',p,'\MRFz');
    cd(path);
       f1 = char(append('n_syN_T1_', strcat(p), '_brain.nii'));
       f2 = char(append('n_syN_T2_', strcat(p), '_brain.nii'));
       f3 = char(append('n_syN_GM_', strcat(p), '_brain.nii'));
       f4 = char(append('n_syN_WM_', strcat(p), '_brain.nii'));
       a = niftiread(f1);
       b = niftiread(f2);
       c = niftiread(f3);
       d = niftiread(f4);
       oma = mT1;
       omb = mT2;
       omc = mGM;
       omd = mWM;
       mT1 = (mT1*(n-1) + a) ./ n;
       mT2 = (mT2*(n-1) + b) ./ n;
       mGM = (mGM*(n-1) + c) ./ n;
       mWM = (mWM*(n-1) + d) ./ n;
       if n > 1
%             sdT1 = sqrt(((n-2)*sdT1.^2+(n-1)/n*(mT1.^2 + a.^2 - 2.*mT1.*a))./(n-1));
%             sdT2 = sqrt(((n-2)*sdT2.^2+(n-1)/n*(mT2.^2 + b.^2 - 2.*mT2.*b))./(n-1));
%             sdGM = sqrt(((n-2)*sdGM.^2+(n-1)/n*(mGM.^2 + c.^2 - 2.*mGM.*c))./(n-1));
%             sdWM = sqrt(((n-2)*sdWM.^2+(n-1)/n*(mWM.^2 + d.^2 - 2.*mWM.*d))./(n-1));           
            sdT1 = sqrt( ((n-2)*sdT1.^2  +  (a-mT1).*(a-oma) )./(n-1));
            sdT2 = sqrt( ((n-2)*sdT2.^2  +  (b-mT2).*(b-omb) )./(n-1));
            sdGM = sqrt( ((n-2)*sdGM.^2  +  (c-mGM).*(c-omc) )./(n-1));
            sdWM = sqrt( ((n-2)*sdWM.^2  +  (d-mWM).*(d-omd) )./(n-1));
       end     
n = n + 1
end


%% get smoothed SD image
    copyfile MNI_T1.nii T1sd.nii;
    copyfile MNI_T2.nii T2sd.nii;
    copyfile MNI_GM_prob.nii GMsd.nii;
    copyfile MNI_WM_prob.nii WMsd.nii;

    as = load_untouch_nii('T1sd.nii');
    bs = load_untouch_nii('T2sd.nii');
    cs = load_untouch_nii('GMsd.nii');
    ds = load_untouch_nii('WMsd.nii');
    
    as.img = sdT1;
    bs.img = sdT2;
    cs.img = sdGM;
    ds.img = sdWM;

    save_untouch_nii(as,'T1sd.nii');
    save_untouch_nii(bs,'T2sd.nii');
    save_untouch_nii(cs,'GMsd.nii');
    save_untouch_nii(ds,'WMsd.nii');

    spm('Defaults','fMRI');
    spm_jobman('initcfg');
    matlabbatch{1}.spm.spatial.smooth.data = {
                                                  char(strcat('T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRFz\T1sd.nii'))
                                                  char(strcat('T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRFz\T2sd.nii'))
                                                  char(strcat('T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRFz\GMsd.nii'))
                                                  char(strcat('T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRFz\WMsd.nii'))
                                                  };
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        
    spm_jobman('run', matlabbatch, cell(0, 8));

as = load_untouch_nii('sT1sd.nii');
bs = load_untouch_nii('sT2sd.nii');
cs = load_untouch_nii('sGMsd.nii');
ds = load_untouch_nii('sWMsd.nii');

sdT1_sm = as.img;
sdT2_sm = bs.img;
sdGM_sm = cs.img;
sdWM_sm = ds.img;

%% generate z score maps
for p = subjID
    path = strcat(MRF_path,'\',p,'\MRFz');
    cd(path);
%     copyfile MNI_T1.nii T1z.nii;
%     copyfile MNI_T2.nii T2z.nii;
%     copyfile MNI_GM_prob.nii GMz.nii;
%     copyfile MNI_WM_prob.nii WMz.nii;
%     copyfile MNI_T1.nii T1z_sm.nii;
%     copyfile MNI_T2.nii T2z_sm.nii;
%     copyfile MNI_GM_prob.nii GMz_sm.nii;
%     copyfile MNI_WM_prob.nii WMz_sm.nii;

% normal SD
    a = load_untouch_nii('T1z.nii');
    b = load_untouch_nii('T2z.nii');
    c = load_untouch_nii('GMz.nii');
    d = load_untouch_nii('WMz.nii');

    a.img = (a.img - mT1)./sdT1;
    b.img = (b.img - mT2)./sdT2;
    c.img = (c.img - mGM)./sdGM;
    d.img = (d.img - mWM)./sdWM;
    
    save_untouch_nii(a,'T1z.nii');
    save_untouch_nii(b,'T2z.nii');
    save_untouch_nii(c,'GMz.nii');
    save_untouch_nii(d,'WMz.nii');

% smoothed SD
    as = load_untouch_nii('T1z_sm.nii');
    bs = load_untouch_nii('T2z_sm.nii');
    cs = load_untouch_nii('GMz_sm.nii');
    ds = load_untouch_nii('WMz_sm.nii');
    
    as.img = (as.img - mT1)./sdT1_sm;
    bs.img = (bs.img - mT2)./sdT2_sm;
    cs.img = (cs.img - mGM)./sdGM_sm;
    ds.img = (ds.img - mWM)./sdWM_sm;
    
    save_untouch_nii(as,'T1z_sm.nii');
    save_untouch_nii(bs,'T2z_sm.nii');
    save_untouch_nii(cs,'GMz_sm.nii');
    save_untouch_nii(ds,'WMz_sm.nii');

end

%% Patient MRFz images
subjID = [XXXX];

MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path);
    copyfile MNI_T1.nii T1z.nii;
    copyfile MNI_T2.nii T2z.nii;
    copyfile MNI_GM_prob.nii GMz.nii;
    copyfile MNI_WM_prob.nii WMz.nii;
    copyfile MNI_T1.nii T1z_sm.nii;
    copyfile MNI_T2.nii T2z_sm.nii;
    copyfile MNI_GM_prob.nii GMz_sm.nii;
    copyfile MNI_WM_prob.nii WMz_sm.nii;

% normal SD
    a = load_untouch_nii('T1z.nii');
    b = load_untouch_nii('T2z.nii');
    c = load_untouch_nii('GMz.nii');
    d = load_untouch_nii('WMz.nii');

    a.img = (a.img - mT1)./sdT1;
    b.img = (b.img - mT2)./sdT2;
    c.img = (c.img - mGM)./sdGM;
    d.img = (d.img - mWM)./sdWM;
    
    save_untouch_nii(a,'T1z.nii');
    save_untouch_nii(b,'T2z.nii');
    save_untouch_nii(c,'GMz.nii');
    save_untouch_nii(d,'WMz.nii');

% smoothed SD
    as = load_untouch_nii('T1z_sm.nii');
    bs = load_untouch_nii('T2z_sm.nii');
    cs = load_untouch_nii('GMz_sm.nii');
    ds = load_untouch_nii('WMz_sm.nii');
    
    as.img = (as.img - mT1)./sdT1_sm;
    bs.img = (bs.img - mT2)./sdT2_sm;
    cs.img = (cs.img - mGM)./sdGM_sm;
    ds.img = (ds.img - mWM)./sdWM_sm;
    
    save_untouch_nii(as,'T1z_sm.nii');
    save_untouch_nii(bs,'T2z_sm.nii');
    save_untouch_nii(cs,'GMz_sm.nii');
    save_untouch_nii(ds,'WMz_sm.nii');

    p

end
        
%-----------------------------------------------------------------------
% Job saved on 18-Nov-2021 16:44:09 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
clear all
% subjID=[XXX];
subjID=[XXX];

spm('Defaults','PET');
spm_jobman('initcfg');
for p = subjID
    p
    matlabbatch{1}.spm.spatial.smooth.data = {
                                              char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\', p, '\MRF_VBM\MNI_T1.nii'))
                                              char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\', p, '\MRF_VBM\MNI_T2.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\', p, '\MRF_VBM\MNI_CSF_fn.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\', p, '\MRF_VBM\MNI_GM_fn.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\', p, '\MRF_VBM\MNI_WM_fn.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\', p, '\MRF_VBM\MNI_extension_fn.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\', p, '\MRF_VBM\MNI_junction_fn.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\', p, '\MRF_VBM\MNI_thickness_fn.nii'))
                                              };
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    
    spm_jobman('run', matlabbatch, cell(0, 2));
end


%% volunteer
subjID=[XX];

subjID=[XXX];
spm('Defaults','PET');
spm_jobman('initcfg');
for p = subjID

    matlabbatch{1}.spm.spatial.smooth.data = {
                                              char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRF_VBM\MNI_T1.nii'))
                                              char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRF_VBM\MNI_T2.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRF_VBM\MNI_CSF_fn.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRF_VBM\MNI_GM_fn.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRF_VBM\MNI_WM_fn.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRF_VBM\MNI_extension_fn.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRF_VBM\MNI_junction_fn.nii'))
%                                               char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRF_VBM\MNI_thickness_fn.nii'))
                                              };
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    
    spm_jobman('run', matlabbatch, cell(0, 8));
end

%% 
clear all
subjID=[XXX];


spm('Defaults','fMRI');
spm_jobman('initcfg');
for p = subjID

    matlabbatch{1}.spm.spatial.smooth.data = {
                                              char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\', p, '\MRF_VBM\MNI_T1z.nii'))
                                              char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\', p, '\MRF_VBM\MNI_T2z.nii'))
                                              };
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    
    spm_jobman('run', matlabbatch, cell(0, 8));
end

%%
% clear all

spm('Defaults','PET');
spm_jobman('initcfg');
matlabbatch{1}.spm.spatial.smooth.data = {
    char(strcat('T:\Imaging\Multimodal\MRF\Peter\MRFzmaps\junction_std.nii'))
    char(strcat('T:\Imaging\Multimodal\MRF\Peter\MRFzmaps\junction_T1_std.nii'))
    char(strcat('T:\Imaging\Multimodal\MRF\Peter\MRFzmaps\junction_T2_std.nii'))
    };
matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

spm_jobman('run', matlabbatch, cell(0, 8));
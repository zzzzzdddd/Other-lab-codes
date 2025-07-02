subjID=[XXXXX];
mT1w = zeros(182,218,182);
sdT1w = zeros(182,218,182);
mT1 = zeros(182,218,182);
sdT1 = zeros(182,218,182);
mT2 = zeros(182,218,182);
sdT2 = zeros(182,218,182);
n_matrix = zeros(182,218,182);
% n_matrix_p = zeros(182,218,182);
% n_matrix_p2 = zeros(182,218,182);
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';

n = 0
for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    p
    f1 = 'MNI_MRFbinary.nii';
    f2 = 'MNI_MRFbinary_T1.nii';
    f3 = 'MNI_MRFbinary_T2.nii';
%     f4 = 'MNI_motion_mask.nii';
    f5 = 'MNI_GM_prob.nii';
    f6 = 'MNI_WM_prob.nii';

%     a = niftiread(f1);
%     b = niftiread(f2);
%     c = niftiread(f3);

    a = double(load_untouch_nii(f1).img);
    b = double(load_untouch_nii(f2).img);
    c = double(load_untouch_nii(f3).img);
    d = zeros(182,218,182);
    d(a>0) = 1;
    e = double(load_untouch_nii(f5).img);
    f = double(load_untouch_nii(f6).img);
    e(e>=0.5) = 1;
    e(e<0.5) = 0;
    f(f>=0.5) = 1;
    f(f<0.5) = 0;
%     imshow(d(:,:,100))
%     d = niftiread(f4);
%     d(d==0&a==0) = 1;
    
%     n_matrix_p2 = n_matrix_p;
    n_matrix_p = n_matrix;
    n_matrix = n_matrix+d.*e;
%     n_p = n;    
%     n = n+1;

%     oma = mT1w;
%     omb = mT1;
%     omc = mT2;

    mT1w = mT1w + a.*d.*e;
    mT1 = mT1 + b.*d.*e;
    mT2 = mT2 + c.*d.*e;

end

mT1w = mT1w./n_matrix;
mT1 = mT1./n_matrix;
mT2 = mT2./n_matrix;
%
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps

file = load_untouch_nii('junction_mean.nii');
file.img = zeros(182,218,182);
file.img = mT1w;
save_untouch_nii(file, 'junction_mean.nii');


file = load_untouch_nii('junction_T1_mean.nii');
file.img = zeros(182,218,182);
file.img = mT1;
save_untouch_nii(file, 'junction_T1_mean.nii');


file = load_untouch_nii('junction_T2_mean.nii');
file.img = zeros(182,218,182);
file.img = mT2;
save_untouch_nii(file, 'junction_T2_mean.nii');

%%
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps

file = load_untouch_nii('T1_mean.nii');
file.img = zeros(182,218,182);
file.img = mT1;
save_untouch_nii(file, 'T1_mean.nii');

file = load_untouch_nii('T2_mean.nii');
file.img = zeros(182,218,182);
file.img = mT2;
save_untouch_nii(file, 'T2_mean.nii');

file = load_untouch_nii('T1_std.nii');
file.img = zeros(182,218,182);
file.img = sdT1;
save_untouch_nii(file, 'T1_std.nii');

file = load_untouch_nii('T2_std.nii');
file.img = zeros(182,218,182);
file.img = sdT2;
save_untouch_nii(file, 'T2_std.nii');
%% smooth the mean map
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps
t1max = 1759; 
t2max = 106;
t1min = 743;
t2min = 29;
T1i = double(load_untouch_nii('T1_mean.nii').img);
T2i = double(load_untouch_nii('T2_mean.nii').img);
kernelSize = 5;
%     [filteredImage1, filteredImage2] = smooth3Dnonzero(T1i, T2i, kernelSize, csfi);
sigma = 6/(2*sqrt(2*log(2)));

% Create a 3D Gaussian kernel
[x, y, z] = meshgrid(-floor(kernelSize/2):floor(kernelSize/2));
gaussianKernel = exp(-(x.^2 + y.^2 + z.^2)/(2*sigma^2));
gaussianKernel = gaussianKernel / sum(gaussianKernel(:));
filteredImage = zeros(size(T1i));

% Apply Gaussian filtering
for i = 1:size(T1i, 1)
    for j = 1:size(T1i, 2)
        for k = 1:size(T1i, 3)
            % Define the neighborhood
            rowRange = max(1, i-floor(kernelSize/2)):min(size(T1i, 1), i+floor(kernelSize/2));
            colRange = max(1, j-floor(kernelSize/2)):min(size(T1i, 2), j+floor(kernelSize/2));
            depthRange = max(1, k-floor(kernelSize/2)):min(size(T1i, 3), k+floor(kernelSize/2));

            % Extract the neighborhood
            neighborhood1 = T1i(rowRange, colRange, depthRange);
            neighborhood2 = T2i(rowRange, colRange, depthRange);

            mask1 = (neighborhood1 < t1max)&(neighborhood1 > t1min);
            mask2 = (neighborhood2 < t2max)&(neighborhood2 > t2min);
            filteredNeighborhood1 = neighborhood1 .* mask1;
            filteredNeighborhood2 = neighborhood2 .* mask2;

            gaussianKernelAdjusted = gaussianKernel(1:length(rowRange), 1:length(colRange), 1:length(depthRange));
            gaussianKernelAdjusted1 = gaussianKernelAdjusted .* mask1;
            gaussianKernelAdjusted2 = gaussianKernelAdjusted .* mask2;
            kernelSum1 = sum(gaussianKernelAdjusted1(:));
            kernelSum2 = sum(gaussianKernelAdjusted2(:));

            if kernelSum1 > 0
                gaussianKernelNormalized1 = gaussianKernelAdjusted1 / kernelSum1;
                filteredValue1 = sum(sum(sum(filteredNeighborhood1 .* gaussianKernelNormalized1)));
            else
                filteredValue1 = T1i(i, j, k);
            end

            if kernelSum2 > 0
                gaussianKernelNormalized2 = gaussianKernelAdjusted2 / kernelSum2;
                filteredValue2 = sum(sum(sum(filteredNeighborhood2 .* gaussianKernelNormalized2)));
            else
                filteredValue2 = T2i(i, j, k);
            end

            % Store the filtered value
            filteredImage1(i, j, k) = filteredValue1;
            filteredImage2(i, j, k) = filteredValue2;

        end
    end
end

newim1 = load_untouch_nii('smT1_nocsf.nii');
newim1.img = filteredImage1;
save_untouch_nii(newim1, 'smT1_nocsf.nii');

newim2 = load_untouch_nii('smT2_nocsf.nii');
newim2.img = filteredImage2;
save_untouch_nii(newim2, 'smT2_nocsf.nii');

%% simple smooth
copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\T1_mean.nii C:\MRFzmaps\T1_mean.nii
copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\T1_std.nii C:\MRFzmaps\T1_std.nii
copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\T2_mean.nii C:\MRFzmaps\T2_mean.nii
copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\T2_std.nii C:\MRFzmaps\T2_std.nii
copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\GM_mean.nii C:\MRFzmaps\GM_mean.nii
copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\GM_std.nii C:\MRFzmaps\GM_std.nii
copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\WM_mean.nii C:\MRFzmaps\WM_mean.nii
copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\WM_std.nii C:\MRFzmaps\WM_std.nii
copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\CSF_mean.nii C:\MRFzmaps\CSF_mean.nii
copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\CSF_std.nii C:\MRFzmaps\CSF_std.nii

spm('Defaults','PET');
spm_jobman('initcfg');
matlabbatch{1}.spm.spatial.smooth.data = {
    char(strcat('C:\MRFzmaps\T1_mean.nii'))
    char(strcat('C:\MRFzmaps\T2_mean.nii'))
    char(strcat('C:\MRFzmaps\T1_std.nii'))
    char(strcat('C:\MRFzmaps\T2_std.nii'))
    char(strcat('C:\MRFzmaps\GM_mean.nii'))
    char(strcat('C:\MRFzmaps\GM_std.nii'))
    char(strcat('C:\MRFzmaps\WM_mean.nii'))
    char(strcat('C:\MRFzmaps\WM_std.nii'))
    char(strcat('C:\MRFzmaps\CSF_mean.nii'))
    char(strcat('C:\MRFzmaps\CSF_std.nii'))
    };
matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

spm_jobman('run', matlabbatch, cell(0, 10));

copyfile C:\MRFzmaps\sT1_mean.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sT1_mean.nii
copyfile C:\MRFzmaps\sT1_std.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sT1_std.nii 
copyfile C:\MRFzmaps\sT2_mean.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sT2_mean.nii
copyfile C:\MRFzmaps\sT2_std.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sT2_std.nii
copyfile C:\MRFzmaps\sGM_mean.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sGM_mean.nii
copyfile C:\MRFzmaps\sGM_std.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sGM_std.nii
copyfile C:\MRFzmaps\sWM_mean.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sWM_mean.nii
copyfile C:\MRFzmaps\sWM_std.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sWM_std.nii
copyfile C:\MRFzmaps\sCSF_mean.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sCSF_mean.nii
copyfile C:\MRFzmaps\sCSF_std.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sCSF_std.nii
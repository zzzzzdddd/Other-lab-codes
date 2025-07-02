copyfile T1.nii T1_255.nii
a = load_untouch_nii('T1_255.nii');
x = single(a.img);
mini = min(x,[],[1 2 3]);
maxi = max(x,[],[1 2 3]);
norm = (x-mini)./(maxi-mini)*255;
a.img = norm;
save_untouch_nii(a,'T1_255.nii')
%%
iterNum_outer=15;  % outer iteration
iterCM=2;  % inner interation for C and M
iter_b=1;  % inner iteration for bias
q = 1.5;   % fuzzifier
th_bg = 5;  %% threshold for removing background
N_region = 3; %% number of tissue types, e.g. WM, GM, CSF
tissueLabel=[1, 2, 3];
str_vector{1} = 'T1_255.nii';   % input a sequence of image file names
% str_vector{2} = 'brainweb_byte_B3N1.nii';
MICO_3Dseq(str_vector, N_region, q, th_bg, iterNum_outer, iter_b, iterCM, tissueLabel);
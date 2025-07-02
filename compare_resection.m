subjID=[XXXX];
for p=subjID
    p
    prefile = char(strcat('Z:\Imaging\Multimodal\Myelin\Patients\', p, '\ROI_T1w.nii'));
    postfile = strcat('C:\Users\eegrvw\Downloads\Resection_analysis\', p, '\Post\Results');
    cd(postfile)
    postfile = dir('Subj_*_ExtractedPostImg_1.nii').name;
    pre = load_untouch_nii(prefile);
    post = load_untouch_nii(postfile);
    A = double(pre.img);
    A = A/max(A,[],'all');
    A(A>0) = 1;
    B = double(post.img);
    B = B/max(B,[],'all');
    B(B>0) = 1;
    tot = sum(A,[1 2 3]);
    p = B - 2*A;
    resect_common = sum(p(p==-1)*(-1));
    percent = resect_common/tot
end
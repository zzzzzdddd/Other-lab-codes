subjID=[XXXXX];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
cd('T:\Imaging\Multimodal\MRF\Peter');
dict = readtable('atlas.xlsx');
atlas = load_untouch_nii('Talairach-labels-1mm.nii');
atlasi = single(atlas.img);
lobes = [];
for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    roi = load_untouch_nii('MNI_ROI.nii');
    roii = single(roi.img);
    roii(roii>0) = 1;
    if p == XXXX
        roii(:,269:end,:) = 0;
    end
%  By center
    [roiblob, nBlobs] = bwlabeln(roii);
    centroid= regionprops3(roiblob, "Centroid");
    c = round(centroid.Centroid);
    cx = c(2);
    cy = c(1);
    cz = c(3);
    labnum = atlasi(cx,cy,cz);

% % by most coverage
%     coms = roii(roii>0).*atlasi(roii>0);
%     labnum = mode(coms);
    if labnum > 0
        brainregion = dict(labnum,2);
        parsed = split(string(brainregion.Var2),".");
        if length(parsed)>1
            lobes = vertcat(lobes, parsed(2));
        else
            lobes = vertcat(lobes, parsed(1));
        end
    end
    
end
[C, ia, ic] = unique(lobes);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts]
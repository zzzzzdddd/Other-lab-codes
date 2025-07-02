%% binarize ANN output
subjID=["P01_12251" "P15_12610" "P16_12669" "P19_12748" "P21_12791" "P22_12836"...
    "P23_12911" "P24_13157" "P31_13375" "P42_13702" "P50_13890" "P62_14212"...
    "P63_14218" "P72_14364" "P76_14439" "P81_14473" "P83_14516" "P85_14705"...
    "P52_13923" "P56_14105" "P57_14129"];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni = single(ann.img);
    anni(anni>0) = 1;
    ann.img = anni;
    save_untouch_nii(ann, 'MNI_lesionprob_ANN.nii');
end


%%
subjID=["P41_13696" "P42_13702" "P50_13890" "P63_14218" "P72_14364" "P76_14439" "P81_14473" ...
    "P83_14516" "P24_13157" "P01_12251" "P02_12253" "P12_12515" "P17_12709" "P43_13706"...
    "P62_14212" "P45_13724" "P52_13923" "P56_14105" "P57_14129" "P68_14287" "P82_14494" "P10_12479"  "P85_14705" "P16_12669"...
    "P15_12610" "P18_12713" "P19_12748" "P21_12791" "P22_12836" "P23_12911" "P32_13395" "P33_13509" "P31_13375" "P61_14203" "P50_13890"];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
t1tp = [];
t2tp = [];
gmtp = [];
wmtp = [];

t1fp = []; 
t2fp = [];
gmfp = [];
wmfp = [];


for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    roi = load_untouch_nii('MNI_ROI_final.nii');
    roii = single(roi.img);
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni =  single(ann.img);
    anni(anni>0) = 1;
    t1 = load_untouch_nii('MNI_T1.nii');
    t1i =  single(t1.img);
    t2 = load_untouch_nii('MNI_T2.nii');
    t2i =  single(t2.img);
    gm = load_untouch_nii('MNI_GM_prob.nii');
    gmi =  single(gm.img);
    wm = load_untouch_nii('MNI_WM_prob.nii');
    wmi =  single(wm.img);
    
    fp = anni-roii;
    fp(fp<0) = 0;
    tp = roii.*anni;

    a = t1i.*tp;
    a = a(a>0);
    t1tp = vertcat(t1tp,a);
    a = t1i.*fp;
    a = a(a>0);
    t1fp = vertcat(t1fp,a);

    a = t2i.*tp;
    a = a(a>0);
    t2tp = vertcat(t2tp,a);
    a = t2i.*fp;
    a = a(a>0);
    t2fp = vertcat(t2fp,a);

    a = gmi.*tp;
    a = a(a>0);
    gmtp = vertcat(gmtp,a);
    a = gmi.*fp;
    a = a(a>0);
    gmfp = vertcat(gmfp,a);

    a = wmi.*tp;
    a = a(a>0);
    wmtp = vertcat(wmtp,a);
    a = wmi.*fp;
    a = a(a>0);
    wmfp = vertcat(wmfp,a);
end

%% compare histogram
[h1,p1] = ttest2(t1tp, t1fp);
[h2,p2] = ttest2(t2tp, t2fp);
[h3,p3] = ttest2(gmtp, gmfp);
[h4,p4] = ttest2(wmtp, wmfp);

figure()
subplot(2,1,1)
histogram(t1tp, 'BinWidth', 5)
title('Distribution of T1 in true positive clusters');
subtitle(strcat("mean = ",string(round(mean(t1tp),2))," std = ",string(round(std(t1tp),2))));
subplot(2,1,2)
histogram(t1fp, 'BinWidth', 5)
title('Distribution of T1 in false positive clusters');
subtitle(strcat("mean = ",string(round(mean(t1fp),2))," std = ",string(round(std(t1fp),2))," p = ",string(p1)));

figure()
subplot(2,1,1)
histogram(t2tp, 'BinWidth', 5)
title('Distribution of T2 in true positive clusters')
subtitle(strcat("mean = ",string(round(mean(t2tp),2))," std = ",string(round(std(t2tp),2))));
xlim([0 200])
subplot(2,1,2)
histogram(t2fp, 'BinWidth', 5)
title('Distribution of T2 in false positive clusters')
subtitle(strcat("mean = ",string(round(mean(t2fp),2))," std = ",string(round(std(t2fp),2))," p = ",string(p2)));
xlim([0 200])

figure()
subplot(2,1,1)
histogram(gmtp, 'BinWidth', 0.05)
title('Distribution of GM prob. in true positive clusters');
subtitle(strcat("mean = ",string(round(mean(gmtp),2))," std = ",string(round(std(gmtp),2))));
subplot(2,1,2)
histogram(gmfp, 'BinWidth', 0.05)
title('Distribution of GM prob. in false positive clusters');
subtitle(strcat("mean = ",string(round(mean(gmfp),2))," std = ",string(round(std(gmfp),2))," p = ",string(p3)));

figure()
subplot(2,1,1)
histogram(wmtp, 'BinWidth', 0.05)
title('Distribution of WM prob. in true positive clusters');
subtitle(strcat("mean = ",string(round(mean(wmtp),2))," std = ",string(round(std(wmtp),2))));
subplot(2,1,2)
histogram(wmfp, 'BinWidth', 0.05)
title('Distribution of WM prob. in false positive clusters');
subtitle(strcat("mean = ",string(round(mean(wmfp),2))," std = ",string(round(std(wmfp),2))," p = ",string(p4)));
%% artifact probability map
subjID=["V01_12506" "V03_12576" "V04_12578" "V05_12579" "V06_12598" "V07_12599" "V08_12605" "V09_12676"...
    "V10_12681" "V11_13486" "V12_13502" "V13_13521" "V14_13522" "V15_13523" "V16_13524" "V17_13528" "V18_13530"...
    "V19_13531" "V20_13534" "V21_13535" "V22_13536" "V23_13538" "V24_13539" "V25_13540" "V26_13543" "V27_13544"...
    "V28_13545" "V29_13555" "V30_13561" "V31_13563" "V32_13575" "V33_13583" "V34_13587" "V35_13591" ...
    "V37_13596" "V38_13603" "V39_13622" "V40_13716" "V41_13738" "V42_13749" "V43_13754" "V44_14271" "V45_14315"...
    "V46_14645" "V47_14649" "V48_14678" "V49_14681"];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni =  single(ann.img);
    if p == "V01_12506"
        artprob = anni;
    else
        artprob = artprob + anni;
    end
end
artprob = artprob./size(subjID,2);
med = median(artprob(artprob>0));
artprob(artprob<med) = 0;
artprob(artprob>0) = 1;

%% generate result after filtering with artifact probability map

subjID=["P41_13696" "P42_13702" "P50_13890" "P63_14218" "P72_14364" "P76_14439" "P81_14473" ...
    "P83_14516" "P24_13157" "P01_12251" "P02_12253" "P12_12515" "P17_12709" "P43_13706"...
    "P62_14212" "P45_13724" "P52_13923" "P56_14105" "P57_14129" "P68_14287" "P82_14494" "P10_12479"  "P85_14705" "P16_12669"...
    "P15_12610" "P18_12713" "P19_12748" "P21_12791" "P22_12836" "P23_12911" "P32_13395" "P33_13509" "P31_13375" "P61_14203" "P50_13890"];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
t1tp = [];
t2tp = [];
gmtp = [];
wmtp = [];

t1fp = []; 
t2fp = [];
gmfp = [];
wmfp = [];


for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    roi = load_untouch_nii('MNI_ROI_final.nii');
    roii = single(roi.img);
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni =  single(ann.img);
    anni(anni>0) = 1;
    anni = anni - artprob;
    anni(anni<0) = 0;
    t1 = load_untouch_nii('MNI_T1.nii');
    t1i =  single(t1.img);
    t2 = load_untouch_nii('MNI_T2.nii');
    t2i =  single(t2.img);
    gm = load_untouch_nii('MNI_GM_prob.nii');
    gmi =  single(gm.img);
    wm = load_untouch_nii('MNI_WM_prob.nii');
    wmi =  single(wm.img);
    
    fp = anni-roii;
    fp(fp<0) = 0;
    tp = roii.*anni;

    a = t1i.*tp;
    a = a(a>0);
    t1tp = vertcat(t1tp,a);
    a = t1i.*fp;
    a = a(a>0);
    t1fp = vertcat(t1fp,a);

    a = t2i.*tp;
    a = a(a>0);
    t2tp = vertcat(t2tp,a);
    a = t2i.*fp;
    a = a(a>0);
    t2fp = vertcat(t2fp,a);

    a = gmi.*tp;
    a = a(a>0);
    gmtp = vertcat(gmtp,a);
    a = gmi.*fp;
    a = a(a>0);
    gmfp = vertcat(gmfp,a);

    a = wmi.*tp;
    a = a(a>0);
    wmtp = vertcat(wmtp,a);
    a = wmi.*fp;
    a = a(a>0);
    wmfp = vertcat(wmfp,a);
end

%% compare histogram
[h1,p1] = ttest2(t1tp, t1fp);
[h2,p2] = ttest2(t2tp, t2fp);
[h3,p3] = ttest2(gmtp, gmfp);
[h4,p4] = ttest2(wmtp, wmfp);

figure()
subplot(2,1,1)
histogram(t1tp, 'BinWidth', 5)
title('Distribution of T1 in true positive clusters');
subtitle(strcat("mean = ",string(round(mean(t1tp),2))," std = ",string(round(std(t1tp),2))));
subplot(2,1,2)
histogram(t1fp, 'BinWidth', 5)
title('Distribution of T1 in false positive clusters');
subtitle(strcat("mean = ",string(round(mean(t1fp),2))," std = ",string(round(std(t1fp),2))," p = ",string(p1)));

figure()
subplot(2,1,1)
histogram(t2tp, 'BinWidth', 5)
title('Distribution of T2 in true positive clusters')
subtitle(strcat("mean = ",string(round(mean(t2tp),2))," std = ",string(round(std(t2tp),2))));
xlim([0 200])
subplot(2,1,2)
histogram(t2fp, 'BinWidth', 5)
title('Distribution of T2 in false positive clusters')
subtitle(strcat("mean = ",string(round(mean(t2fp),2))," std = ",string(round(std(t2fp),2))," p = ",string(p2)));
xlim([0 200])

figure()
subplot(2,1,1)
histogram(gmtp, 'BinWidth', 0.05)
title('Distribution of GM prob. in true positive clusters');
subtitle(strcat("mean = ",string(round(mean(gmtp),2))," std = ",string(round(std(gmtp),2))));
subplot(2,1,2)
histogram(gmfp, 'BinWidth', 0.05)
title('Distribution of GM prob. in false positive clusters');
subtitle(strcat("mean = ",string(round(mean(gmfp),2))," std = ",string(round(std(gmfp),2))," p = ",string(p3)));

figure()
subplot(2,1,1)
histogram(wmtp, 'BinWidth', 0.05)
title('Distribution of WM prob. in true positive clusters');
subtitle(strcat("mean = ",string(round(mean(wmtp),2))," std = ",string(round(std(wmtp),2))));
subplot(2,1,2)
histogram(wmfp, 'BinWidth', 0.05)
title('Distribution of WM prob. in false positive clusters');
subtitle(strcat("mean = ",string(round(mean(wmfp),2))," std = ",string(round(std(wmfp),2))," p = ",string(p4)));

%% blob extraction 
% function binaryImage = ExtractNLargestBlobs(binaryImage, numberToExtract)
% 	[labeledImage, numberOfBlobs] = bwlabeln(binaryImage);
% 	blobMeasurements = regionprops3(labeledImage, 'area');
% 	allAreas = [blobMeasurements.Area];
% 	if numberToExtract > 0
% 		[sortedAreas, sortIndexes] = sort(allAreas, 'descend');
% 	end
% 	% Extract the "numberToExtract" largest blob(a)s using ismember().
% 	if  numberOfBlobs < numberToExtract
%         biggestBlob = ismember(labeledImage, sortIndexes(1:numberOfBlobs));
%     else
%         biggestBlob = ismember(labeledImage, sortIndexes(1:numberToExtract));
%     end
% 	% Convert from integer labeled image into binary (logical) image.
% 	binaryImage = biggestBlob > 0;
% end

subjID=["P01_12251" "P15_12610" "P16_12669" "P19_12748" "P21_12791" "P22_12836"...
    "P23_12911" "P24_13157" "P31_13375" "P42_13702" "P50_13890" "P62_14212"...
    "P63_14218" "P72_14364" "P76_14439" "P81_14473" "P83_14516" "P85_14705"...
    "P52_13923" "P56_14105" "P57_14129"];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

xtall = [];
ytall = [];
xfall = [];
yfall = [];
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni = single(ann.img);
    [labeledImage, numberOfBlobs] = bwlabeln(anni);
	blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
	allinds = [blobMeasurements.VoxelIdxList];

    roi = load_untouch_nii('MNI_ROI_final.nii');
    roii = single(roi.img);
    [roiblob, nBlobs] = bwlabeln(roii);
    roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
    roiinds = roiMeasurements.VoxelIdxList{1};
    
    trueind = [];
    for i = 1:numberOfBlobs
        A = blobMeasurements.VoxelIdxList{i};
        if sum(ismember(roiinds,A))>0
            trueind = horzcat(trueind, i);
        end
    end

    t1s = load_untouch_nii('MNI_T1.nii');
    t1si = single(t1s.img);
    t2s = load_untouch_nii('MNI_T2.nii');
    t2si = single(t2s.img);
    

    xt = [];
    yt = [];
    xf = [];
    yf = [];
    for i = 1:numberOfBlobs
        if ismember(i,trueind) == 1
            xt = vertcat(xt, t1si(allinds{i}));
            yt = vertcat(yt, t2si(allinds{i}));
     
        else
            xf = vertcat(xf, t1si(allinds{i}));
            yf = vertcat(yf, t2si(allinds{i}));
        end
    end

    figure()
    scatter(xt, yt, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
    hold on
    scatter(xf, yf, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
    title(p)
    xlabel('T1(ms)');
    ylabel('T2(ms)');

    xtall = vertcat(xtall, xt);
    ytall = vertcat(ytall, yt);
    xfall = vertcat(xfall, xf);
    yfall = vertcat(yfall, yf);
end

figure()
scatter(xtall, ytall, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
hold on
scatter(xfall, yfall, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
title('all patients')
xlabel('T1(ms)');
ylabel('T2(ms)');

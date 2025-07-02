%% get mean MRF T1 and T2 from controls, skip if already has data
subjID=["V01_12506" "V03_12576" "V04_12578" "V05_12579" "V06_12598" "V07_12599" "V08_12605" "V09_12676"...
    "V10_12681" "V11_13486" "V12_13502" "V13_13521" "V14_13522" "V15_13523" "V16_13524" "V17_13528" "V18_13530"...
    "V19_13531" "V20_13534" "V21_13535" "V22_13536" "V23_13538" "V24_13539" "V25_13540" "V26_13543" "V27_13544"...
    "V28_13545" "V29_13555" "V30_13561" "V31_13563" "V32_13575" "V33_13583" "V34_13587" "V35_13591" ...
    "V37_13596" "V38_13603" "V39_13622" "V40_13716" "V41_13738" "V42_13749" "V43_13754" "V44_14271" "V45_14315"...
    "V46_14645" "V47_14649" "V48_14678" "V49_14681"];

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
%% Promising patients and separate into GM and WM
subjID=["P15_12610" "P16_12669" "P19_12748" "P21_12791" "P22_12836"...
    "P23_12911" "P24_13157" "P31_13375" "P62_14212"...
    "P83_14516" "P85_14705"...
    "P52_13923" "P56_14105" "P57_14129"];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

xtallW = [];
ytallW = [];
xfallW = [];
yfallW = [];
xtallG = [];
ytallG = [];
xfallG = [];
yfallG = [];
t1max = 2000; 
t2max = 180;
t1min = 600;
t2min = 20;
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni = single(ann.img);
    anni(anni>0) = 1;
    mask = single(load_untouch_nii('MNI_Brain_Mask.nii').img);
    mask(mask>0.95) = 1;
    anni = anni.*mask;
    GM = single(load_untouch_nii('MNI_GM.nii').img);
    GM(GM<0.6) = 0;
    GM(GM>=0.6) = 1;
    WM = single(load_untouch_nii('MNI_WM.nii').img);
    WM(WM<0.6) = 0;
    WM(WM>=0.6) = 1;
    
    % GM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni.*GM);
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
    

    xtGM = [];
    ytGM = [];
    xfGM = [];
    yfGM = [];
    for i = 1:numberOfBlobs
        if ismember(i,trueind) == 1
            xtGM = vertcat(xtGM, t1si(allinds{i}));
            ytGM = vertcat(ytGM, t2si(allinds{i}));
        else
            xfGM = vertcat(xfGM, t1si(allinds{i}));
            yfGM = vertcat(yfGM, t2si(allinds{i}));
        end
    end
    xtGMc = xtGM;
    ytGMc = ytGM;
    xfGMc = xfGM;
    yfGMc = yfGM;
    xtGM((xtGMc > t1max)|(xtGMc< t1min)|(ytGMc > t2max)|(ytGMc< t2min)) = []; 
    ytGM((ytGMc > t2max)|(ytGMc< t2min)|(xtGMc > t1max)|(xtGMc< t1min)) = [];
    xfGM((xfGMc > t1max)|(xfGMc< t1min)|(yfGMc > t2max)|(yfGMc< t2min)) = []; 
    yfGM((yfGMc > t2max)|(yfGMc< t2min)|(xfGMc > t1max)|(xfGMc< t1min)) = [];


    figure()
    scatter(xtGM, ytGM, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
    hold on
    scatter(xfGM, yfGM, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
    title(p,"GM")
    xlabel('T1(ms)');
    ylabel('T2(ms)');

    % WM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni.*WM);
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
    

    xtWM = [];
    ytWM = [];
    xfWM = [];
    yfWM = [];

    for i = 1:numberOfBlobs
        if ismember(i,trueind) == 1
            xtWM = vertcat(xtWM, t1si(allinds{i}));
            ytWM = vertcat(ytWM, t2si(allinds{i}));
        else
            xfWM = vertcat(xfWM, t1si(allinds{i}));
            yfWM = vertcat(yfWM, t2si(allinds{i}));
        end
    end
    xtWMc = xtWM;
    ytWMc = ytWM;
    xfWMc = xfWM;
    yfWMc = yfWM;
    xtWM((xtWMc > t1max)|(xtWMc< t1min)|(ytWMc > t2max)|(ytWMc< t2min)) = []; 
    ytWM((ytWMc > t2max)|(ytWMc< t2min)|(xtWMc > t1max)|(xtWMc< t1min)) = [];
    xfWM((xfWMc > t1max)|(xfWMc< t1min)|(yfWMc > t2max)|(yfWMc< t2min)) = []; 
    yfWM((yfWMc > t2max)|(yfWMc< t2min)|(xfWMc > t1max)|(xfWMc< t1min)) = [];

    figure()
    scatter(xtWM, ytWM, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
    hold on
    scatter(xfWM, yfWM, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
    title(p,"WM")
    xlabel('T1(ms)');
    ylabel('T2(ms)');

    xtallW = vertcat(xtallW, xtWM);
    ytallW = vertcat(ytallW, ytWM);
    xfallW = vertcat(xfallW, xfWM);
    yfallW = vertcat(yfallW, yfWM);

    xtallG = vertcat(xtallG, xtGM);
    ytallG = vertcat(ytallG, ytGM);
    xfallG = vertcat(xfallG, xfGM);
    yfallG = vertcat(yfallG, yfGM);
end

% figure()
% scatter(xtallW, ytallW, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% hold on
% scatter(xfallW, yfallW, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% title('all patients WM')
% xlabel('T1(ms)');
% ylabel('T2(ms)');
% xlim([0 1800]);
% ylim([0 100]);
% 
% figure()
% scatter(xtallG, ytallG, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% hold on
% scatter(xfallG, yfallG, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% title('all patients GM')
% xlabel('T1(ms)');
% ylabel('T2(ms)');
% xlim([0 1800]);
% ylim([0 100]);

% calculate mean and std 
%white matter
meanwt1t = round(mean(xtallW));
stdwt1t = round(std(xtallW));
meanwt2t = round(mean(ytallW),1);
stdwt2t = round(std(ytallW),1);
meanwt1f = round(mean(xfallW));
stdwt1f = round(std(xfallW));
meanwt2f = round(mean(yfallW),1);
stdwt2f = round(std(yfallW),1);
[h1,p1w] = ttest2(xtallW, xfallW,'Vartype','unequal');
[h2,p2w] = ttest2(ytallW, yfallW,'Vartype','unequal');
%gray matter
meangt1t = round(mean(xtallG));
stdgt1t = round(std(xtallG));
meangt2t = round(mean(ytallG),1);
stdgt2t = round(std(ytallG),1);
meangt1f = round(mean(xfallG));
stdgt1f = round(std(xfallG));
meangt2f = round(mean(yfallG),1);
stdgt2f = round(std(yfallG),1);
[h1,p1g] = ttest2(xtallG, xfallG,'Vartype','unequal');
[h2,p2g] = ttest2(ytallG, yfallG,'Vartype','unequal');
% make tables
GM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meangt1t),string(char(177)),string(stdgt1t)); strcat(string(meangt1f),string(char(177)),string(stdgt1f));string(p1g)];
T2 = [strcat(string(meangt2t),string(char(177)),string(stdgt2t)); strcat(string(meangt2f),string(char(177)),string(stdgt2f));string(p2g)];
GMtab1 = table(GM,T1,T2)

WM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meanwt1t),string(char(177)),string(stdwt1t)); strcat(string(meanwt1f),string(char(177)),string(stdwt1f));string(p1w)];
T2 = [strcat(string(meanwt2t),string(char(177)),string(stdwt2t)); strcat(string(meanwt2f),string(char(177)),string(stdwt2f));string(p2w)];
WMtab1 = table(WM,T1,T2)
%% Promising patients and separate into GM and WM probability maps non-normalized
subjID=["P15_12610" "P16_12669" "P19_12748" "P21_12791" "P22_12836"...
    "P23_12911" "P24_13157" "P31_13375" "P62_14212"...
    "P83_14516" "P85_14705"...
    "P52_13923" "P56_14105" "P57_14129"];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

xtallW = [];
ytallW = [];
xfallW = [];
yfallW = [];
xtallG = [];
ytallG = [];
xfallG = [];
yfallG = [];
t1max = 2000; 
t2max = 180;
t1min = 600;
t2min = 20;
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni = single(ann.img);
    anni(anni>0) = 1;
    mask = single(load_untouch_nii('MNI_Brain_Mask.nii').img);
    mask(mask>0.95) = 1;
    anni = anni.*mask;
    GM = single(load_untouch_nii('MNI_GM_prob.nii').img);
    GM(GM<0.8) = 0;
    GM(GM>=0.8) = 1;
    WM = single(load_untouch_nii('MNI_WM_prob.nii').img);
    WM(WM<0.8) = 0;
    WM(WM>=0.8) = 1;
    
    % GM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni.*GM);
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
    

    xtGM = [];
    ytGM = [];
    xfGM = [];
    yfGM = [];
    for i = 1:numberOfBlobs
        if ismember(i,trueind) == 1
            xtGM = vertcat(xtGM, t1si(allinds{i}));
            ytGM = vertcat(ytGM, t2si(allinds{i}));
        else
            xfGM = vertcat(xfGM, t1si(allinds{i}));
            yfGM = vertcat(yfGM, t2si(allinds{i}));
        end
    end
    xtGMc = xtGM;
    ytGMc = ytGM;
    xfGMc = xfGM;
    yfGMc = yfGM;
    xtGM((xtGMc > t1max)|(xtGMc< t1min)|(ytGMc > t2max)|(ytGMc< t2min)) = []; 
    ytGM((ytGMc > t2max)|(ytGMc< t2min)|(xtGMc > t1max)|(xtGMc< t1min)) = [];
    xfGM((xfGMc > t1max)|(xfGMc< t1min)|(yfGMc > t2max)|(yfGMc< t2min)) = []; 
    yfGM((yfGMc > t2max)|(yfGMc< t2min)|(xfGMc > t1max)|(xfGMc< t1min)) = [];


%     figure()
%     scatter(xtGM, ytGM, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     hold on
%     scatter(xfGM, yfGM, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     title(p,"GM")
%     xlabel('T1(ms)');
%     ylabel('T2(ms)');

    % WM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni.*WM);
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
    

    xtWM = [];
    ytWM = [];
    xfWM = [];
    yfWM = [];

    for i = 1:numberOfBlobs
        if ismember(i,trueind) == 1
            xtWM = vertcat(xtWM, t1si(allinds{i}));
            ytWM = vertcat(ytWM, t2si(allinds{i}));
        else
            xfWM = vertcat(xfWM, t1si(allinds{i}));
            yfWM = vertcat(yfWM, t2si(allinds{i}));
        end
    end
    xtWMc = xtWM;
    ytWMc = ytWM;
    xfWMc = xfWM;
    yfWMc = yfWM;
    xtWM((xtWMc > t1max)|(xtWMc< t1min)|(ytWMc > t2max)|(ytWMc< t2min)) = []; 
    ytWM((ytWMc > t2max)|(ytWMc< t2min)|(xtWMc > t1max)|(xtWMc< t1min)) = [];
    xfWM((xfWMc > t1max)|(xfWMc< t1min)|(yfWMc > t2max)|(yfWMc< t2min)) = []; 
    yfWM((yfWMc > t2max)|(yfWMc< t2min)|(xfWMc > t1max)|(xfWMc< t1min)) = [];

%     figure()
%     scatter(xtWM, ytWM, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     hold on
%     scatter(xfWM, yfWM, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     title(p,"WM")
%     xlabel('T1(ms)');
%     ylabel('T2(ms)');

    xtallW = vertcat(xtallW, xtWM);
    ytallW = vertcat(ytallW, ytWM);
    xfallW = vertcat(xfallW, xfWM);
    yfallW = vertcat(yfallW, yfWM);

    xtallG = vertcat(xtallG, xtGM);
    ytallG = vertcat(ytallG, ytGM);
    xfallG = vertcat(xfallG, xfGM);
    yfallG = vertcat(yfallG, yfGM);
end

% figure()
% scatter(xtallW, ytallW, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% hold on
% scatter(xfallW, yfallW, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% title('all patients WM')
% xlabel('T1(ms)');
% ylabel('T2(ms)');
% xlim([0 2000]);
% ylim([0 400]);

% figure()
% scatter(xtallG, ytallG, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% hold on
% scatter(xfallG, yfallG, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% title('all patients GM')
% xlabel('T1(ms)');
% ylabel('T2(ms)');
% xlim([0 2000]);
% ylim([0 400]);
%% calculate mean and std 
%white matter
meanwt1t = round(mean(xtallW));
stdwt1t = round(std(xtallW));
meanwt2t = round(mean(ytallW),1);
stdwt2t = round(std(ytallW),1);
meanwt1f = round(mean(xfallW));
stdwt1f = round(std(xfallW));
meanwt2f = round(mean(yfallW),1);
stdwt2f = round(std(yfallW),1);
[h1,p1w] = ttest2(xtallW, xfallW,'Vartype','unequal');
[h2,p2w] = ttest2(ytallW, yfallW,'Vartype','unequal');
%gray matter
meangt1t = round(mean(xtallG));
stdgt1t = round(std(xtallG));
meangt2t = round(mean(ytallG),1);
stdgt2t = round(std(ytallG),1);
meangt1f = round(mean(xfallG));
stdgt1f = round(std(xfallG));
meangt2f = round(mean(yfallG),1);
stdgt2f = round(std(yfallG),1);
[h1,p1g] = ttest2(xtallG, xfallG,'Vartype','unequal');
[h2,p2g] = ttest2(ytallG, yfallG,'Vartype','unequal');
% make tables
GM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meangt1t),string(char(177)),string(stdgt1t)); strcat(string(meangt1f),string(char(177)),string(stdgt1f));string(p1g)];
T2 = [strcat(string(meangt2t),string(char(177)),string(stdgt2t)); strcat(string(meangt2f),string(char(177)),string(stdgt2f));string(p2g)];
GMtab2 = table(GM,T1,T2)

WM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meanwt1t),string(char(177)),string(stdwt1t)); strcat(string(meanwt1f),string(char(177)),string(stdwt1f));string(p1w)];
T2 = [strcat(string(meanwt2t),string(char(177)),string(stdwt2t)); strcat(string(meanwt2f),string(char(177)),string(stdwt2f));string(p2w)];
WMtab2 = table(WM,T1,T2)



%% Promising patients and separate into GM and WM normalized
subjID=["P15_12610" "P16_12669" "P19_12748" "P21_12791" "P22_12836"...
    "P23_12911" "P24_13157" "P31_13375" "P62_14212"...
    "P83_14516" "P85_14705"...
    "P52_13923" "P56_14105" "P57_14129"];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

xtallW = [];
ytallW = [];
xfallW = [];
yfallW = [];
xtallG = [];
ytallG = [];
xfallG = [];
yfallG = [];
t1max = 2000; 
t2max = 180;
t1min = 600;
t2min = 20;
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni = single(ann.img);
    mask = single(load_untouch_nii('MNI_Brain_Mask.nii').img);
    anni(anni>0) = 1;
    mask(mask>0.95) = 1;
    anni = anni.*mask;
    GM = single(load_untouch_nii('MNI_GM.nii').img);
    GM(GM<0.8) = 0;
    GM(GM>=0.8) = 1;
    WM = single(load_untouch_nii('MNI_WM.nii').img);
    WM(WM<0.8) = 0;
    WM(WM>=0.8) = 1;
    
    % GM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni.*GM);
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
    

    xtGM = [];
    ytGM = [];
    xfGM = [];
    yfGM = [];
    xtGMn = [];
    ytGMn = [];
    xfGMn = [];
    yfGMn = [];
    for i = 1:numberOfBlobs
        if ismember(i,trueind) == 1
            xtGM = vertcat(xtGM, t1si(allinds{i}));
            ytGM = vertcat(ytGM, t2si(allinds{i}));
            xtGMn = vertcat(xtGMn, t1si(allinds{i})./mT1(allinds{i}));
            ytGMn = vertcat(ytGMn, t2si(allinds{i})./mT2(allinds{i}));
        else
            xfGM = vertcat(xfGM, t1si(allinds{i}));
            yfGM = vertcat(yfGM, t2si(allinds{i}));
            xfGMn = vertcat(xfGMn, t1si(allinds{i})./mT1(allinds{i}));
            yfGMn = vertcat(yfGMn, t2si(allinds{i})./mT2(allinds{i}));
        end
    end
    xtGMc = xtGM;
    ytGMc = ytGM;
    xfGMc = xfGM;
    yfGMc = yfGM;
    xtGMn((xtGMc > t1max)|(xtGMc< t1min)|(ytGMc > t2max)|(ytGMc< t2min)) = []; 
    ytGMn((ytGMc > t2max)|(ytGMc< t2min)|(xtGMc > t1max)|(xtGMc< t1min)) = [];
    xfGMn((xfGMc > t1max)|(xfGMc< t1min)|(yfGMc > t2max)|(yfGMc< t2min)) = []; 
    yfGMn((yfGMc > t2max)|(yfGMc< t2min)|(xfGMc > t1max)|(xfGMc< t1min)) = [];


%     figure()
%     scatter(xtGM, ytGM, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     hold on
%     scatter(xfGM, yfGM, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     title(p,"GM")
%     xlabel('T1(ms)');
%     ylabel('T2(ms)');

    % WM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni.*WM);
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
    

    xtWM = [];
    ytWM = [];
    xfWM = [];
    yfWM = [];
    xtWMn = [];
    ytWMn = [];
    xfWMn = [];
    yfWMn = [];
    for i = 1:numberOfBlobs
        if ismember(i,trueind) == 1
            xtWM = vertcat(xtWM, t1si(allinds{i}));
            ytWM = vertcat(ytWM, t2si(allinds{i}));
            xtWMn = vertcat(xtWMn, t1si(allinds{i})./mT1(allinds{i}));
            ytWMn = vertcat(ytWMn, t2si(allinds{i})./mT2(allinds{i}));
        else
            xfWM = vertcat(xfWM, t1si(allinds{i}));
            yfWM = vertcat(yfWM, t2si(allinds{i}));
            xfWMn = vertcat(xfWMn, t1si(allinds{i})./mT1(allinds{i}));
            yfWMn = vertcat(yfWMn, t2si(allinds{i})./mT2(allinds{i}));
        end
    end
    xtWMc = xtWM;
    ytWMc = ytWM;
    xfWMc = xfWM;
    yfWMc = yfWM;
    xtWMn((xtWMc > t1max)|(xtWMc< t1min)|(ytWMc > t2max)|(ytWMc< t2min)) = []; 
    ytWMn((ytWMc > t2max)|(ytWMc< t2min)|(xtWMc > t1max)|(xtWMc< t1min)) = [];
    xfWMn((xfWMc > t1max)|(xfWMc< t1min)|(yfWMc > t2max)|(yfWMc< t2min)) = []; 
    yfWMn((yfWMc > t2max)|(yfWMc< t2min)|(xfWMc > t1max)|(xfWMc< t1min)) = [];

    xtallW = vertcat(xtallW, xtWMn);
    ytallW = vertcat(ytallW, ytWMn);
    xfallW = vertcat(xfallW, xfWMn);
    yfallW = vertcat(yfallW, yfWMn);

    xtallG = vertcat(xtallG, xtGMn);
    ytallG = vertcat(ytallG, ytGMn);
    xfallG = vertcat(xfallG, xfGMn);
    yfallG = vertcat(yfallG, yfGMn);
end
figure()
scatter(xtallW, ytallW, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
hold on
scatter(xfallW, yfallW, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
title("WM")
xlabel('T1(ms)');
ylabel('T2(ms)');

figure()
scatter(xtallG, ytallG, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
hold on
scatter(xfallG, yfallG, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
title("GM")
xlabel('T1(ms)');
ylabel('T2(ms)');
%% calculate mean and std 
%white matter
meanwt1t = round(mean(xtallW),3);
stdwt1t = round(std(xtallW),3);
meanwt2t = round(mean(ytallW),3);
stdwt2t = round(std(ytallW),3);
meanwt1f = round(mean(xfallW),3);
stdwt1f = round(std(xfallW),3);
meanwt2f = round(mean(yfallW),3);
stdwt2f = round(std(yfallW),3);
[h1,p1w] = ttest2(xtallW, xfallW,'Vartype','unequal');
[h2,p2w] = ttest2(ytallW, yfallW,'Vartype','unequal');
%gray matter
meangt1t = round(mean(xtallG),3);
stdgt1t = round(std(xtallG),3);
meangt2t = round(mean(ytallG),3);
stdgt2t = round(std(ytallG),3);
meangt1f = round(mean(xfallG),3);
stdgt1f = round(std(xfallG),3);
meangt2f = round(mean(yfallG),3);
stdgt2f = round(std(yfallG),3);
[h1,p1g] = ttest2(xtallG, xfallG,'Vartype','unequal');
[h2,p2g] = ttest2(ytallG, yfallG,'Vartype','unequal');
% make tables
GM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meangt1t),string(char(177)),string(stdgt1t)); strcat(string(meangt1f),string(char(177)),string(stdgt1f));string(p1g)];
T2 = [strcat(string(meangt2t),string(char(177)),string(stdgt2t)); strcat(string(meangt2f),string(char(177)),string(stdgt2f));string(p2g)];
GMtab3 = table(GM,T1,T2)

WM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meanwt1t),string(char(177)),string(stdwt1t)); strcat(string(meanwt1f),string(char(177)),string(stdwt1f));string(p1w)];
T2 = [strcat(string(meanwt2t),string(char(177)),string(stdwt2t)); strcat(string(meanwt2f),string(char(177)),string(stdwt2f));string(p2w)];
WMtab3 = table(WM,T1,T2)

%% Promising patients and separate into GM_prob and WM_prob normalized
subjID=["P15_12610" "P16_12669" "P19_12748" "P21_12791" "P22_12836"...
    "P23_12911" "P24_13157" "P31_13375" "P62_14212"...
    "P83_14516" "P85_14705"...
    "P52_13923" "P56_14105" "P57_14129"];
MRF_path='T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

xtallW = [];
ytallW = [];
xfallW = [];
yfallW = [];
xtallG = [];
ytallG = [];
xfallG = [];
yfallG = [];
t1max = 2000; 
t2max = 180;
t1min = 600;
t2min = 20;
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni = single(ann.img);
    anni(anni>0) = 1;
    mask = single(load_untouch_nii('MNI_Brain_Mask.nii').img);
    mask(mask>0.95) = 1;
    anni = anni.*mask;
    GM = single(load_untouch_nii('MNI_GM.nii').img);
    GM(GM<0.8) = 0;
    GM(GM>=0.8) = 1;
    WM = single(load_untouch_nii('MNI_WM.nii').img);
    WM(WM<0.8) = 0;
    WM(WM>=0.8) = 1;
    
    % GM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni.*GM);
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
    

    xtGM = [];
    ytGM = [];
    xfGM = [];
    yfGM = [];
    xtGMn = [];
    ytGMn = [];
    xfGMn = [];
    yfGMn = [];
    for i = 1:numberOfBlobs
        if ismember(i,trueind) == 1
            xtGM = vertcat(xtGM, t1si(allinds{i}));
            ytGM = vertcat(ytGM, t2si(allinds{i}));
            xtGMn = vertcat(xtGMn, t1si(allinds{i})./mT1(allinds{i}));
            ytGMn = vertcat(ytGMn, t2si(allinds{i})./mT2(allinds{i}));
        else
            xfGM = vertcat(xfGM, t1si(allinds{i}));
            yfGM = vertcat(yfGM, t2si(allinds{i}));
            xfGMn = vertcat(xfGMn, t1si(allinds{i})./mT1(allinds{i}));
            yfGMn = vertcat(yfGMn, t2si(allinds{i})./mT2(allinds{i}));
        end
    end
    xtGMc = xtGM;
    ytGMc = ytGM;
    xfGMc = xfGM;
    yfGMc = yfGM;
    xtGMn((xtGMc > t1max)|(xtGMc< t1min)|(ytGMc > t2max)|(ytGMc< t2min)) = []; 
    ytGMn((ytGMc > t2max)|(ytGMc< t2min)|(xtGMc > t1max)|(xtGMc< t1min)) = [];
    xfGMn((xfGMc > t1max)|(xfGMc< t1min)|(yfGMc > t2max)|(yfGMc< t2min)) = []; 
    yfGMn((yfGMc > t2max)|(yfGMc< t2min)|(xfGMc > t1max)|(xfGMc< t1min)) = [];


%     figure()
%     scatter(xtGM, ytGM, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     hold on
%     scatter(xfGM, yfGM, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     title(p,"GM")
%     xlabel('T1(ms)');
%     ylabel('T2(ms)');

    % WM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni.*WM);
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
    

    xtWM = [];
    ytWM = [];
    xfWM = [];
    yfWM = [];
    xtWMn = [];
    ytWMn = [];
    xfWMn = [];
    yfWMn = [];
    for i = 1:numberOfBlobs
        if ismember(i,trueind) == 1
            xtWM = vertcat(xtWM, t1si(allinds{i}));
            ytWM = vertcat(ytWM, t2si(allinds{i}));
            xtWMn = vertcat(xtWMn, t1si(allinds{i})./mT1(allinds{i}));
            ytWMn = vertcat(ytWMn, t2si(allinds{i})./mT2(allinds{i}));
        else
            xfWM = vertcat(xfWM, t1si(allinds{i}));
            yfWM = vertcat(yfWM, t2si(allinds{i}));
            xfWMn = vertcat(xfWMn, t1si(allinds{i})./mT1(allinds{i}));
            yfWMn = vertcat(yfWMn, t2si(allinds{i})./mT2(allinds{i}));
        end
    end
    xtWMc = xtWM;
    ytWMc = ytWM;
    xfWMc = xfWM;
    yfWMc = yfWM;
    xtWMn((xtWMc > t1max)|(xtWMc< t1min)|(ytWMc > t2max)|(ytWMc< t2min)) = []; 
    ytWMn((ytWMc > t2max)|(ytWMc< t2min)|(xtWMc > t1max)|(xtWMc< t1min)) = [];
    xfWMn((xfWMc > t1max)|(xfWMc< t1min)|(yfWMc > t2max)|(yfWMc< t2min)) = []; 
    yfWMn((yfWMc > t2max)|(yfWMc< t2min)|(xfWMc > t1max)|(xfWMc< t1min)) = [];

%     figure()
%     scatter(xtWM, ytWM, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     hold on
%     scatter(xfWM, yfWM, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     title(p,"WM")
%     xlabel('T1(ms)');
%     ylabel('T2(ms)');

    xtallW = vertcat(xtallW, xtWMn);
    ytallW = vertcat(ytallW, ytWMn);
    xfallW = vertcat(xfallW, xfWMn);
    yfallW = vertcat(yfallW, yfWMn);

    xtallG = vertcat(xtallG, xtGMn);
    ytallG = vertcat(ytallG, ytGMn);
    xfallG = vertcat(xfallG, xfGMn);
    yfallG = vertcat(yfallG, yfGMn);
end
figure()
scatter(xtallW, ytallW, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
hold on
scatter(xfallW, yfallW, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
title("WM")
xlabel('T1(ms)');
ylabel('T2(ms)');

figure()
scatter(xtallG, ytallG, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
hold on
scatter(xfallG, yfallG, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
title("GM")
xlabel('T1(ms)');
ylabel('T2(ms)');
% calculate mean and std 
%white matter
meanwt1t = round(mean(xtallW),3);
stdwt1t = round(std(xtallW),3);
meanwt2t = round(mean(ytallW),3);
stdwt2t = round(std(ytallW),3);
meanwt1f = round(mean(xfallW),3);
stdwt1f = round(std(xfallW),3);
meanwt2f = round(mean(yfallW),3);
stdwt2f = round(std(yfallW),3);
[h1,p1w] = ttest2(xtallW, xfallW,'Vartype','unequal');
[h2,p2w] = ttest2(ytallW, yfallW,'Vartype','unequal');
%gray matter
meangt1t = round(mean(xtallG),3);
stdgt1t = round(std(xtallG),3);
meangt2t = round(mean(ytallG),3);
stdgt2t = round(std(ytallG),3);
meangt1f = round(mean(xfallG),3);
stdgt1f = round(std(xfallG),3);
meangt2f = round(mean(yfallG),3);
stdgt2f = round(std(yfallG),3);
[h1,p1g] = ttest2(xtallG, xfallG,'Vartype','unequal');
[h2,p2g] = ttest2(ytallG, yfallG,'Vartype','unequal');
% make tables
GM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meangt1t),string(char(177)),string(stdgt1t)); strcat(string(meangt1f),string(char(177)),string(stdgt1f));string(p1g)];
T2 = [strcat(string(meangt2t),string(char(177)),string(stdgt2t)); strcat(string(meangt2f),string(char(177)),string(stdgt2f));string(p2g)];
GMtab4 = table(GM,T1,T2)

WM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meanwt1t),string(char(177)),string(stdwt1t)); strcat(string(meanwt1f),string(char(177)),string(stdwt1f));string(p1w)];
T2 = [strcat(string(meanwt2t),string(char(177)),string(stdwt2t)); strcat(string(meanwt2f),string(char(177)),string(stdwt2f));string(p2w)];
WMtab4 = table(WM,T1,T2)
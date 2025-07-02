%% previously s117  P89_14878
%% now use s126 P57_14129
cd T:\Imaging\Multimodal\MRF\Peter
load T1cm.mat
load T2cm.mat
s = 126;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\Regis_files_wo_Q_bet
SYNf = load_untouch_nii('XXX.nii');
syn = single(SYNf.img);
syn = syn(:,:,s)';
syn = flip(syn,1);

cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM
MPf = load_untouch_nii('MNI_T1w.nii');
MP = single(MPf.img);
MP = MP(:,:,s)';
MP = flip(MP,1);

roif = load_untouch_nii('MNI_ROI_final.nii');
roi = single(roif.img);
roi = roi(:,:,s)';
roi = flip(roi,1);
[B, L] = bwboundaries(roi,'noholes');

annf = load_untouch_nii('MNI_lesionprob_ANN.nii');
ann = single(annf.img);
ann = ann(:,:,s)';
ann = flip(ann,1);

GMf = load_untouch_nii('MNI_T1z.nii');
GM = single(GMf.img);
GM = GM(:,:,s)';
GM = flip(GM,1);

WMf = load_untouch_nii('MNI_T2z.nii');
WM = single(WMf.img);
WM = WM(:,:,s)';
WM = flip(WM,1);

T1f = load_untouch_nii('MNI_T1.nii');
T1 = single(T1f.img);
T1 = T1(:,:,s)';
T1 = flip(T1,1);

T2f = load_untouch_nii('MNI_T2.nii');
T2 = single(T2f.img);
T2 = T2(:,:,s)';
T2 = flip(T2,1);

predf = load_untouch_nii('prediction.nii');
pred = single(predf.img);
pred = T2(:,:,s)';
pred = flip(T2,1);

annf = load_untouch_nii('MNI_lesionprob_ANN.nii');
ann = single(annf.img);
annfp = ann(:,:,sfp)';
annfp = flip(annfp,1);
anntp = ann(:,:,stp)';
anntp = flip(anntp,1);

roi_manual = T1./T1;


%% Overlay 2 different colors

cd T:\Imaging\Multimodal\MRF\Peter
% sfp = 111;
% stp = 132;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\P57_14129\MRF_VBM


sfp = 82;
stp = 82;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM

sfp = 90;
stp = 90;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM

sfp = 42;
stp = 42;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM

sfp = 112;
stp = 112;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM

% sfp = 68;
% stp = 68;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM

% sfp = 100;
% stp = 100;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM

% sfp = 90;
% stp = 90;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM

% sfp = 111;
% stp = 111;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM
% 
% sfp = 126;
% stp = 126;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM
cmap = [0 0 0
    0.85 0.225 0.1
    0.85 0.25 0.1
    0.85 0.275 0.1
    0.85 0.3 0.1
    0.85 0.325 0.1
    0.85 0.35 0.1
    0.85 0.375 0.1
    0.9 0.4 0.1
    0.9 0.425 0.1
    0.9 0.45 0.1
    0.9 0.475 0.1
    0.9 0.5 0.1
    0.9 0.525 0.1
    0.9 0.55 0.1
    0.9 0.575 0.1
    0.9 0.6 0.1
    0.9 0.625 0.1
    0.9 0.65 0.1
    0.9 0.675 0.1
    0.9 0.7 0.1]; 

MPf = load_untouch_nii('MNI_T1w.nii');
MP = single(MPf.img)*2000;
MPfp = MP(:,:,sfp)';
MPfp = flip(MPfp,1);
MPtp = MP(:,:,stp)';
MPtp = flip(MPtp,1);

roif = load_untouch_nii('MNI_ROI_final.nii');
roi = single(roif.img);
roifp = roi(:,:,sfp)';
roifp = flip(roifp,1);
roitp = roi(:,:,stp)';
roitp = flip(roitp,1);
comsize = size(MPfp);
comsize(2) = comsize(2)*2+10;
combined = zeros(comsize);
combined(:,193:end) = roitp;
[B, L] = bwboundaries(combined,'noholes');

annf = load_untouch_nii('softmax.nii');
ann = single(annf.img);
ann(ann<0.5) = 0;
annfp = ann(:,:,sfp)';
annfp = flip(annfp,1);
anntp = ann(:,:,stp)';
anntp = flip(anntp,1);

figure()
combined1 = combined;
combined2 = combined;

ax1 = axes;
combined1(:,1:182) = annfp;
combined1(:,183:192) = 1;
combined1(:,193:end) = anntp;
A = imagesc(combined1)
colormap(ax1, cmap)
% colorbar;
A.AlphaData = 1;
hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
end

ax2 = axes;
combined2(:,1:182) = MPfp;
combined2(:,193:end) = MPtp;
C = imagesc(combined2)
colormap(ax2, 'gray')
C.AlphaData = 0.7;


linkaxes([ax1,ax2]);
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];



% load T1cm.mat
% load T2cm.mat
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM
% 
% SYNf = load_untouch_nii('T1w_data_brain.nii');
% syn = single(SYNf.img);
% syn = syn(86:222,77:242,s)';
% % syn = syn(:,:,s)';
% syn = flip(syn,1);
% 
% MPf = load_untouch_nii('T1_cat.nii');
% MP = single(MPf.img);
% MP = MP(90:428,41:426,122)';
% MP = flip(MP,1);
% 
% roif = load_untouch_nii('FCD_ROI_final.nii');
% roi = single(roif.img);
% roi = roi(86:222,77:242,s)';
% roi = flip(roi,1);
% [B, L] = bwboundaries(roi,'noholes');


% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM
% annf = load_untouch_nii('MNI_lesionprob_ANN.nii');
% ann = single(annf.img);
% ann = ann(:,:,s)';
% ann = flip(ann,1);
% 
% GMf = load_untouch_nii('MNI_GM.nii');
% GM = single(GMf.img);
% GM = GM(:,:,s)';
% GM = flip(GM,1);
% 
% WMf = load_untouch_nii('MNI_WM.nii');
% WM = single(WMf.img);
% WM = WM(:,:,s)';
% WM = flip(WM,1);

roi_manual = T1./T1;
%% make figures


figure()
subplot(2,4,1)
imshow(MP/1000)
title('T1w MPRAGE')
hold on
for k = 1:length(B)
   boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
end

subplot(2,4,5)
% figure()
imshow(syn)
title('Synthetic T1w')
hold on
for k = 1:length(B)
   boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
end

subplot(2,4,2)
% figure()
imshow(GM)
title('Gray Matter')
% colorbar
% colormap gray

subplot(2,4,6)
% figure()
imshow(WM)
title('White Matter')
% colorbar
% colormap gray

subplot(2,4,3)
% figure()
imshow(T1.*roi_manual,[0 3000])
colormap(T1colormap)
colorbar
title('MRF T1')
hold on
for k = 1:length(B)
   boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
end

subplot(2,4,7)
% figure()
imshow(T2.*roi_manual,[0 200])
colormap(T2colormap)
colorbar
title('MRF T2')
hold on
for k = 1:length(B)
   boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
end

% figure()
subplot(2,4,4)
imshow(roi)
title('Lesion label')

% figure()
subplot(2,4,8)
imshow(ann*10)
title('FCD probability')
colorbar
colormap gray
% hold on
% for k = 1:length(B)
%    boundary = B{k};
% %    plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
% end

%% tiled layout
t = tiledlayout(2,3,'TileSpacing','Compact');
ax = gobjects(1,6);
for i = 1:6
    ax(i) = nexttile();
    switch i
        case 1
            mprage = imshow(syn/3, [])
            ax1 = mprage;
%             ax1.AlphaData = 0.3
            title('T1w MPRAGE', 'FontSize', 15)           
            hold on

%             for k = 1:length(B)
%                 boundary = B{k};
%                 plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
%             end
%             hold all;
            hb = imshow(ann, [])
            ax2 = hb
            ax2.AlphaData = 0.3
            ax2.Visible = 'off'; 

        case 2
%             imshow(GM)
%             title('Gray Matter', 'FontSize', 15)
%             colormap(ax(i),'gray')
%             colorbar
            imshow(roi)
        case 3
%             imshow(WM)
%             title('White Matter', 'FontSize', 15)
%             colormap(ax(i), 'gray')
%             colorbar
            mprage = imshow(MP/3, [])
            ax1 = mprage;
            title('T1w MPRAGE', 'FontSize', 15)           
            hold on
%         case 4
%             imshow(roi)
%             title('Lesion Label', 'FontSize', 15)
%             imshow(syn)
%             title('Synthetic T1w', 'FontSize', 15)
%             hold on
%             for k = 1:length(B)
%                 boundary = B{k};
%                 plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
%             end
        case 4
            imshow(T1.*roi_manual,[0 3000])
            colormap(ax(i), T1colormap)
            colorbar
            title('MRF T1', 'FontSize', 15)
            hold on
            for k = 1:length(B)
                boundary = B{k};
%                 plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
            end
        case 5
            imshow(T2.*roi_manual,[0 200],'Colormap',T2colormap)
            colormap(ax(i), T2colormap)
            colorbar
            title('MRF T2', 'FontSize', 15)
            hold on
            for k = 1:length(B)
                boundary = B{k};
%                 plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
            end
        case 6
            imshow(ann*1.2)
            title('FCD probability', 'FontSize', 15)
%             colormap(ax(i), 'hot')
            colormap(ax(i), 'gray')
            colorbar
            hold on
%             for k = 1:length(B)
%                 boundary = B{k};
%                 plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
%             end


    end
end

%% stack volunteer MRF
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\PXXXX\MRF_VBM
M1 = zeros(182,218,5);
M2 = zeros(182,218,5);

subjID=["XXXXXX"];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
% MRF_path= '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Normal';

for p = 1:size(subjID,2)
    p
    path = strcat(MRF_path,'\',subjID(p),'\MRF_VBM');
    cd(path)
    z1 = load_untouch_nii('MNI_T1.nii');
    t1slice = single(z1.img);
    M1(:,:,p) = t1slice(:,:,130);

    z2 = load_untouch_nii('MNI_T2.nii');
    t2slice = double(z2.img);
    M2(:,:,p) = t2slice(:,:,130);

    z3 = load_untouch_nii('MNI_GM_prob.nii');
    t3slice = single(z3.img);
    M3(:,:,p) = t3slice(:,:,130);

    z4 = load_untouch_nii('MNI_WM_prob.nii');
    t4slice = double(z4.img);
    M4(:,:,p) = t4slice(:,:,130);

end

t = tiledlayout(1,4);
ax = gobjects(1,2);

ax(1)=nexttile()
hs = slice(M1,[],[],1:5);
title('Healthy Control T1')
colormap(ax(1), T1colormap);
shading interp
set(gca,'visible','off')

ax(2)=nexttile()
hs2 = slice(M2,[],[],1:5);
title('Healthy Control T2')
colormap(ax(2), T2colormap);
shading interp
set(gca,'visible','off')

ax(3)=nexttile()
hs = slice(M3,[],[],1:5);
title('Healthy Control T1')
colormap(ax(3), 'gray');
shading interp
set(gca,'visible','off')

ax(4)=nexttile()
hs = slice(M4,[],[],1:5);
title('Healthy Control T2')
colormap(ax(4), 'gray');
shading interp
set(gca,'visible','off')

%% read stuff and better layout
cd T:\Imaging\Multimodal\MRF\Peter
load T1cm.mat
load T2cm.mat
sfp = 116;
stp = 116;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\PXXXX\MRF_VBM
cd Z:\Imaging\Multimodal\Myelin\Patients\XXX
% MPf = load_untouch_nii('MNI_MPRAGE.nii');
MPf = load_untouch_nii('T1.nii');
MP = single(MPf.img);
MPfp = MP(:,:,sfp)';
MPfp = flip(MPfp,1);
MPtp = MP(:,:,stp)';
MPtp = flip(MPtp,1);

% roif = load_untouch_nii('MNI_ROI_final.nii');
roif = load_untouch_nii('ROI_T1w.nii');
roi = single(roif.img);
roifp = roi(:,:,sfp)';
roifp = flip(roifp,1);
roitp = roi(:,:,stp)';
roitp = flip(roitp,1);
% comsize = size(MPfp);
% comsize(2) = comsize(2)*2+10;
% combined = zeros(comsize);
% combined(:,193:end) = roitp;
[B, L] = bwboundaries(combined,'noholes');

% annf = load_untouch_nii('MNI_lesionprob_ANN.nii');
annf = load_untouch_nii('MNI_cortical_mask.nii');
ann = single(annf.img);
annfp = ann(:,:,sfp)';
annfp = flip(annfp,1);
anntp = ann(:,:,stp)';
anntp = flip(anntp,1);

GMf = load_untouch_nii('T1_brain_pve_1.nii');
GM = single(GMf.img);
GMfp = GM(:,:,sfp)';
GMfp = flip(GMfp,1);
GMtp = GM(:,:,stp)';
GMtp = flip(GMtp,1);

WMf = load_untouch_nii('T1_brain_pve_2.nii');
WM = single(WMf.img);
WMfp = WM(:,:,sfp)';
WMfp = flip(WMfp,1);
WMtp = WM(:,:,stp)';
WMtp = flip(WMtp,1);

T1f = load_untouch_nii('MRF_T1.nii');
T1 = single(T1f.img);
T1fp = T1(:,:,sfp)';
T1fp = flip(T1fp,1);
T1tp = T1(:,:,stp)';
T1tp = flip(T1tp,1);

T2f = load_untouch_nii('MRF_T2.nii');
T2 = single(T2f.img);
T2fp = T2(:,:,sfp)';
T2fp = flip(T2fp,1);
T2tp = T2(:,:,stp)';
T2tp = flip(T2tp,1);

roi_manual = T1./T1;


t = tiledlayout(2,3,'TileSpacing','tight','Padding','tight', ...
    'Position', [0.1300 0.1100 0.1 0.1],'OuterPosition', [0 0 1 1]);
ax = gobjects(1,6);
for i = 1:6
    ax(i) = nexttile();
    comsize = size(MPfp);
    comsize(2) = comsize(2)*2+10;
%     combined = zeros(comsize);
    switch i
        case 1            
%             combined(:,1:182) = MPfp/1000;
%             combined(:,183:192) = 1000;
%             combined(:,193:end) = MPtp/1000;
            combined = MPfp/1000;
            mprage = imshow(combined, [])
            ax1 = mprage;
            ax1.AlphaData = 0.3
            title('T1w MPRAGE', 'FontSize', 15)
            xlabel('Slice 111                                   Slice 132', 'FontSize', 15)

        case 2
%             combined(:,1:182) = GMfp;
%             combined(:,183:192) = 1;
%             combined(:,193:end) = GMtp;
            combined = GMfp;
            imshow(combined)
            title('Gray Matter', 'FontSize', 15)
            colormap(ax(i),'gray')
            colorbar
            xlabel('Slice 111                                   Slice 132', 'FontSize', 15)
        case 3
%             combined(:,1:182) = WMfp;
%             combined(:,183:192) = 1;
%             combined(:,193:end) = WMtp;
            combined = WMfp;
            imshow(combined)
            title('White Matter', 'FontSize', 15)
            colormap(ax(i),'gray')
            colorbar         
            xlabel('Slice 111                                   Slice 132', 'FontSize', 15)
        case 4
%             combined(:,1:182) = T1fp.*roi_manual(:,:,sfp)';
%             combined(:,183:192) = 10000;
%             combined(:,193:end) = T1tp.*roi_manual(:,:,stp)';   
            combined = T1fp;
            imshow(combined,[0 3000])
            colormap(ax(i), T1colormap);
            title('MRF T1', 'FontSize', 15)
            colorbar
            xlabel('Slice 111                                   Slice 132', 'FontSize', 15)
%             hold on
%             for k = 1:length(B)
%                 boundary = B{k};
%                 plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
%             end        


        case 5
%             combined(:,1:182) = T2fp.*roi_manual(:,:,sfp)';
%             combined(:,183:192) = 10000;
%             combined(:,193:end) = T2tp.*roi_manual(:,:,stp)';
            combined = T2fp
            imshow(combined,[0 200])
            colormap(ax(i), T2colormap);
            title('MRF T2', 'FontSize', 15)
            colorbar
            xlabel('Slice 111                                   Slice 132', 'FontSize', 15)
%             hold on
%             for k = 1:length(B)
%                 boundary = B{k};
%                 plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
%             end

        case 6
%             combined(:,1:182) = annfp;
%             combined(:,183:192) = 1;
%             combined(:,193:end) = anntp;
            combined = annfp;
            imshow(combined)
            colormap(ax(i), 'hot');
            title('FCD probability', 'FontSize', 15)
            colorbar
            xlabel('Slice 111                                   Slice 132', 'FontSize', 15)
            hold on
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
            end

    end
end
%% Overlay 2 different colors
figure()
combined1 = combined;
combined2 = combined;

ax1 = axes;
combined1(:,1:182) = annfp;
combined1(:,183:192) = 1;
combined1(:,193:end) = anntp;
A = imagesc(combined1)
colormap(ax1, 'hot')
A.AlphaData = 1;

ax2 = axes;
combined2(:,1:182) = MPfp;
combined2(:,183:192) = 100;
combined2(:,193:end) = MPtp;
C = imagesc(combined2)
colormap(ax2, 'gray')
C.AlphaData = 0.7;

hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
end
% linkaxes([ax1,ax2]);
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% set([ax1, ax2],'color','none','visible','off');
% over = ind2rgb(ann*100, 'hot')

%% 14364 slice 116 for myelin project
cd T:\Imaging\Multimodal\MRF\Peter
load T1cm.mat
load T2cm.mat
s = 116;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\PXXX\Regis_files_wo_Q_bet
% SYNf = load_untouch_nii('n_syN_T1w_PXXXX_Warped.nii');
% syn = single(SYNf.img);
% syn = syn(:,:,s)';
% syn = flip(syn,1);

cd Z:\Imaging\Multimodal\Myelin\Patients\PXXX
MPf = load_untouch_nii('T1.nii');
MP = single(MPf.img);
MP = MP(:,:,s)';
MP = flip(MP,1);

roif = load_untouch_nii('ROI_T1w.nii');
roi = single(roif.img);
roi = roi(:,:,s)';
roi = flip(roi,1);
[B, L] = bwboundaries(roi,'noholes');

annf = load_untouch_nii('T1oT2_nmCC.nii.nii');
ann = single(annf.img);
ann = ann(:,:,s)';
ann = flip(ann,1);

% GMf = load_untouch_nii('MNI_GM.nii');
% GM = single(GMf.img);
% GM = GM(:,:,s)';
% GM = flip(GM,1);
% 
% WMf = load_untouch_nii('MNI_WM.nii');
% WM = single(WMf.img);
% WM = WM(:,:,s)';
% WM = flip(WM,1);

T1f = load_untouch_nii('MRF_T1.nii');
T1 = single(T1f.img);
T1 = T1(:,:,s)';
T1 = flip(T1,1);

T2f = load_untouch_nii('MRF_T2.nii');
T2 = single(T2f.img);
T2 = T2(:,:,s)';
T2 = flip(T2,1);

roi_manual = T1./T1;

%%
t = tiledlayout(1,4,'TileSpacing','tight','Padding','tight', ...
    'Position', [0.1300 0.1100 0.1 0.1],'OuterPosition', [0 0 1 1]);
ax = gobjects(1,4);
for i = 1:4
    ax(i) = nexttile();
    switch i
        case 1            
            mprage = imshow(MP, [])
            ax1 = mprage;
            title('T1w MPRAGE', 'FontSize', 15)
            hold on
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
            end  
        
        case 2
            t1ot2image = imshow(ann.*10000000, [])
            ax2 = t1ot2image;
            title('T1w/T2w', 'FontSize', 15)

        case 3
            imshow(T1,[0 3000])
            colormap(ax(i), T1colormap);
            title('MRF T1', 'FontSize', 15)
            colorbar
      
        case 4
            imshow(T2,[0 200])
            colormap(ax(i), T2colormap);
            title('MRF T2', 'FontSize', 15)
            colorbar
    end
end


cd T:\Imaging\Multimodal\MRF\Peter
sfp = 111;
stp = 132;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\PXXXX\MRF_VBM
MPf = load_untouch_nii('MNI_MPRAGE.nii');
MP = single(MPf.img);
MPfp = MP(:,:,sfp)';
MPfp = flip(MPfp,1);
MPtp = MP(:,:,stp)';
MPtp = flip(MPtp,1);

roif = load_untouch_nii('MNI_ROI_final.nii');
roi = single(roif.img);
roifp = roi(:,:,sfp)';
roifp = flip(roifp,1);
roitp = roi(:,:,stp)';
roitp = flip(roitp,1);
comsize = size(MPfp);
comsize(2) = comsize(2)*2+10;
combined = zeros(comsize);
combined(:,193:end) = roitp;
[B, L] = bwboundaries(combined,'noholes');

annf = load_untouch_nii('prediction.nii');
ann = single(annf.img);
annfp = ann(:,:,sfp)';
annfp = flip(annfp,1);
anntp = ann(:,:,stp)';
anntp = flip(anntp,1);

figure()
combined1 = combined;
combined2 = combined;

ax1 = axes;
combined1(:,1:182) = annfp;
combined1(:,183:192) = 1;
combined1(:,193:end) = anntp;
A = imagesc(combined1)
colormap(ax1, 'warm')
A.AlphaData = 1;

ax2 = axes;
combined2(:,1:182) = MPfp;
combined2(:,183:192) = 100;
combined2(:,193:end) = MPtp;
C = imagesc(combined2)
colormap(ax2, 'gray')
C.AlphaData = 0.7;

hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 1)
end
% linkaxes([ax1,ax2]);
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% set([ax1, ax2],'color','none','visible','off');
% over = ind2rgb(ann*100, 'hot')%% non-mni space 
% previously P89 at 104
% P57 at 96

load T1cm.mat
load T2cm.mat
s = 97;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\PXXXX
SYNf = load_untouch_nii('T1w_data_brain.nii');
syn = single(SYNf.img);
syn = syn(86:222,77:242,s)';
% syn = syn(:,:,s)';
syn = flip(syn,1);

MPf = load_untouch_nii('T1_cat.nii');
MP = single(MPf.img);
MP = MP(90:428,41:426,122)';
MP = flip(MP,1);

roif = load_untouch_nii('FCD_ROI_final.nii');
roi = single(roif.img);
roi = roi(86:222,77:242,s)';
roi = flip(roi,1);
[B, L] = bwboundaries(roi,'noholes');


% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\P57_14129\MRF_VBM
% annf = load_untouch_nii('MNI_lesionprob_ANN.nii');
% ann = single(annf.img);
% ann = ann(:,:,s)';
% ann = flip(ann,1);
% 
% GMf = load_untouch_nii('MNI_GM.nii');
% GM = single(GMf.img);
% GM = GM(:,:,s)';
% GM = flip(GM,1);
% 
% WMf = load_untouch_nii('MNI_WM.nii');
% WM = single(WMf.img);
% WM = WM(:,:,s)';
% WM = flip(WM,1);

roi_manual = T1./T1;

%% nnUNet paper example figure
cd T:\Imaging\Multimodal\MRF\Peter
load T1cm.mat
load T2cm.mat
% sfp = 130;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM
% sfp = 108;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM
%     ax(i).XLim = [20 160];
%     ax(i).YLim = [25 190];
% sfp = 132;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM
%     ax(i).XLim = [28 152];
%     ax(i).YLim = [43 177];
sfp = 35;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM
sfp = 112;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM

% ax(i).XLim = [30 150];
% ax(i).YLim = [50 190];
% sfp = 112;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM
%     ax(i).XLim = [20 160];
%     ax(i).YLim = [30 185];
MPf = load_untouch_nii('MNI_T1w.nii');
MP = single(MPf.img);
MPfp = MP(:,:,sfp)';
MPfp = flip(MPfp,1);

roif = load_untouch_nii('MNI_ROI_final.nii');
roi = single(roif.img);
roifp = roi(:,:,sfp)';
roifp = flip(roifp,1);

combined = roifp;
[B, L] = bwboundaries(combined,'noholes');

t1cf = load_untouch_nii('MNI_MPRAGE_brain.nii');
t1c = single(t1cf.img);
t1cfp = t1c(:,:,sfp)';
t1cfp = flip(t1cfp,1);

GMf = load_untouch_nii('MNI_GM_prob.nii');
GM = single(GMf.img);
GMfp = GM(:,:,sfp)';
GMfp = flip(GMfp,1);

WMf = load_untouch_nii('MNI_WM_prob.nii');
WM = single(WMf.img);
WMfp = WM(:,:,sfp)';
WMfp = flip(WMfp,1);

CSFf = load_untouch_nii('MNI_CSF_prob.nii');
CSF = single(CSFf.img);
CSFfp = CSF(:,:,sfp)';
CSFfp = flip(CSFfp,1);

T1f = load_untouch_nii('MNI_T1.nii');
T1 = single(T1f.img);
T1fp = T1(:,:,sfp)';
T1fp = flip(T1fp,1);

T2f = load_untouch_nii('MNI_T2.nii');
T2 = single(T2f.img);
T2fp = T2(:,:,sfp)';
T2fp = flip(T2fp,1);

roi_manual = T1./T1;

t = tiledlayout(2,4,'TileSpacing','tight','Padding','tight', ...
    'Position', [0.1300 0.1100 0.1 0.1],'OuterPosition', [0 0 1 1]);
ax = gobjects(2,4);
cmap = [0 0 0
    0.85 0.225 0.1
    0.85 0.25 0.1
    0.85 0.275 0.1
    0.85 0.3 0.1
    0.85 0.325 0.1
    0.85 0.35 0.1
    0.85 0.375 0.1
    0.9 0.4 0.1
    0.9 0.425 0.1
    0.9 0.45 0.1
    0.9 0.475 0.1
    0.9 0.5 0.1
    0.9 0.525 0.1
    0.9 0.55 0.1
    0.9 0.575 0.1
    0.9 0.6 0.1
    0.9 0.625 0.1
    0.9 0.65 0.1
    0.9 0.675 0.1
    0.9 0.7 0.1]; 
for i = 1:8
    ax(i) = nexttile();
    comsize = size(MPfp);
    combined = zeros(comsize);
    switch i
        case 1            
            combined(:,1:182) = t1cfp;
            mprage = imshow(combined, [])
            ax1 = mprage;
            title('T1w MPRAGE', 'FontSize', 15)
            hold on
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
            end
        case 2            
            combined(:,1:182) = MPfp*2000;
            mprage = imshow(combined, [])
            ax1 = mprage;
            title('T1w syn', 'FontSize', 15)
            hold on
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
            end
        case 3
            combined(:,1:182) = T1fp;
            imshow(combined,[0 3000])
            colormap(ax(i), T1colormap);
            title('MRF T1', 'FontSize', 15)
%             colorbar
            hold on
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
            end        
        case 4
            combined(:,1:182) = T2fp;
            imshow(combined,[0 200])
            colormap(ax(i), T2colormap);
            title('MRF T2', 'FontSize', 15)
%             colorbar
            hold on
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
            end
        case 5
            combined(:,1:182) = T2fp;
%             imshow(combined,[0 200])
%             colormap(ax(i), T2colormap);
            combined(:,1:182) = T1fp;
            imshow(combined,[0 3000])
            colormap(ax(i), T1colormap);
%             title('MRF T1', 'FontSize', 15)
%             colorbar
            hold on
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
            end
        case 6
            combined(:,1:182) = GMfp;
            imshow(combined,[])
            title('GM', 'FontSize', 15)
            hold on
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
            end         
        case 7
            combined(:,1:182) = WMfp;
            imshow(combined,[])
            title('WM', 'FontSize', 15)
            hold on
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
            end     
        case 8
            combined(:,1:182) = CSFfp;
            imshow(combined,[])
            title('CSF', 'FontSize', 15)
            hold on
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
            end      
    end
ax(i).XLim = [30 150];
ax(i).YLim = [50 190];
end

%% nnUNet paper example figure
cd T:\Imaging\Multimodal\MRF\Peter
load T1cm.mat
load T2cm.mat
% sfp = 130;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM
sfp = 108;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM
%     ax(i).XLim = [20 160];
%     ax(i).YLim = [25 190];
% sfp = 132;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM
%     ax(i).XLim = [28 152];
%     ax(i).YLim = [43 177];
% sfp = 35;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM
% ax(i).XLim = [30 150];
% ax(i).YLim = [50 190];
% sfp = 112;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM
%     ax(i).XLim = [20 160];
%     ax(i).YLim = [30 185];
MPf = load_untouch_nii('MNI_T1w.nii');
MP = single(MPf.img);
MPfp = MP(:,:,sfp)';
MPfp = flip(MPfp,1);

roif = load_untouch_nii('MNI_ROI_final.nii');
roi = single(roif.img);
roifp = roi(:,:,sfp)';
roifp = flip(roifp,1);

combined = roifp;
[B, L] = bwboundaries(combined,'noholes');

t1cf = load_untouch_nii('MNI_MPRAGE_brain.nii');
t1c = single(t1cf.img);
t1cfp = t1c(:,:,sfp)';
t1cfp = flip(t1cfp,1);

GMf = load_untouch_nii('MNI_GM_prob.nii');
GM = single(GMf.img);
GMfp = GM(:,:,sfp)';
GMfp = flip(GMfp,1);

WMf = load_untouch_nii('MNI_WM_prob.nii');
WM = single(WMf.img);
WMfp = WM(:,:,sfp)';
WMfp = flip(WMfp,1);

CSFf = load_untouch_nii('MNI_CSF_prob.nii');
CSF = single(CSFf.img);
CSFfp = CSF(:,:,sfp)';
CSFfp = flip(CSFfp,1);

T1f = load_untouch_nii('MNI_T1.nii');
T1 = single(T1f.img);
T1fp = T1(:,:,sfp)';
T1fp = flip(T1fp,1);

T2f = load_untouch_nii('MNI_T2.nii');
T2 = single(T2f.img);
T2fp = T2(:,:,sfp)';
T2fp = flip(T2fp,1);

roi_manual = T1./T1;

t = tiledlayout(2,4,'TileSpacing','tight','Padding','tight', ...
    'Position', [0.1300 0.1100 0.1 0.1],'OuterPosition', [0 0 1 1]);
ax = gobjects(2,4);
cmap = [0 0 0
    0.85 0.225 0.1
    0.85 0.25 0.1
    0.85 0.275 0.1
    0.85 0.3 0.1
    0.85 0.325 0.1
    0.85 0.35 0.1
    0.85 0.375 0.1
    0.9 0.4 0.1
    0.9 0.425 0.1
    0.9 0.45 0.1
    0.9 0.475 0.1
    0.9 0.5 0.1
    0.9 0.525 0.1
    0.9 0.55 0.1
    0.9 0.575 0.1
    0.9 0.6 0.1
    0.9 0.625 0.1
    0.9 0.65 0.1
    0.9 0.675 0.1
    0.9 0.7 0.1]; 
for i = 1:8
    ax(i) = nexttile();
    comsize = size(MPfp);
    combined = zeros(comsize);
    switch i
        case 1            
            combined(:,1:182) = t1cfp;
            mprage = imshow(combined, [])
            ax1 = mprage;
            title('T1w MPRAGE', 'FontSize', 15)
        case 2            
            combined(:,1:182) = MPfp*2000;
            mprage = imshow(combined, [])
            ax1 = mprage;
            title('T1w syn', 'FontSize', 15)
        case 3
            combined(:,1:182) = T1fp;
            imshow(combined,[0 3000])
            colormap(ax(i), T1colormap);
            title('MRF T1', 'FontSize', 15)
%             colorbar
      
        case 4
            combined(:,1:182) = T2fp;
            imshow(combined,[0 200])
            colormap(ax(i), T2colormap);
            title('MRF T2', 'FontSize', 15)
%             colorbar

        case 5
            combined(:,1:182) = T2fp;
%             imshow(combined,[0 200])
%             colormap(ax(i), T2colormap);
            combined(:,1:182) = T1fp;
            imshow(combined,[0 3000])
            colormap(ax(i), T1colormap);
%             title('MRF T1', 'FontSize', 15)
%             colorbar
        case 6
            combined(:,1:182) = GMfp;
            imshow(combined,[])
            title('GM', 'FontSize', 15)      
        case 7
            combined(:,1:182) = WMfp;
            imshow(combined,[])
            title('WM', 'FontSize', 15)   
        case 8
            combined(:,1:182) = CSFfp;
            imshow(combined,[])
            title('CSF', 'FontSize', 15)    
    end
    ax(i).XLim = [20 160];
    ax(i).YLim = [25 190];
end

%% show case MRF FP vs VBM FP
cd T:\Imaging\Multimodal\MRF\Peter
% sfp = 111;
% stp = 132;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM

sfp = 112;
stp = 110;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM
% 
% sfp = 115;
% stp = 114;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXX\MRF_VBM
% 
% sfp = 80;
% stp = 80;
% cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM

cmap = [0 0 0
    0.85 0.225 0.1
    0.85 0.25 0.1
    0.85 0.275 0.1
    0.85 0.3 0.1
    0.85 0.325 0.1
    0.85 0.35 0.1
    0.85 0.375 0.1
    0.9 0.4 0.1
    0.9 0.425 0.1
    0.9 0.45 0.1
    0.9 0.475 0.1
    0.9 0.5 0.1
    0.9 0.525 0.1
    0.9 0.55 0.1
    0.9 0.575 0.1
    0.9 0.6 0.1
    0.9 0.625 0.1
    0.9 0.65 0.1
    0.9 0.675 0.1
    0.9 0.7 0.1]; 

cmap2 = [0 0 0
    1 0.6 0.85
    1 0.675 0.85
    1 0.65 0.85
    1 0.625 0.85
    1 0.6 0.85
    1 0.575 0.85
    1 0.55 0.85
    1 0.525 0.85
    1 0.5 0.85
    1 0.475 0.85
    1 0.45 0.85
    1 0.425 0.85
    1 0.4 0.85
    1 0.375 0.85
    1 0.35 0.85
    1 0.325 0.85
    1 0.3 0.85
    1 0.275 0.85
    1 0.25 0.85
    1 0.225 0.85]; 


MPf = load_untouch_nii('MNI_T1w.nii');
MP = single(MPf.img)*2000;
MPfp = MP(:,:,sfp)';
MPfp = flip(MPfp,1);
MPtp = MP(:,:,stp)';
MPtp = flip(MPtp,1);

roif = load_untouch_nii('MNI_ROI_final.nii');
roi = single(roif.img);
roifp = roi(:,:,sfp)';
roifp = flip(roifp,1);
roitp = roi(:,:,stp)';
roitp = flip(roitp,1);
comsize = size(MPfp);
comsize(2) = comsize(2)*2+10;
combined = zeros(comsize);
combined(:,193:end) = roitp;
[B, L] = bwboundaries(combined,'noholes');

annf = load_untouch_nii('softmax.nii');
ann = single(annf.img);
ann(ann<0.05) = 0;
annfp = ann(:,:,sfp)';
annfp = flip(annfp,1);
anntp = ann(:,:,stp)';
anntp = flip(anntp,1);

valnnf = load_untouch_nii('MNI_lesionprob_ANN_fn.nii');
valnn = single(valnnf.img);
valnn(valnn<0.05) = 0;
valnnfp = valnn(:,:,sfp)';
valnnfp = flip(valnnfp,1);
valnntp = valnn(:,:,stp)';
valnntp = flip(valnntp,1);

figure()
combined1 = combined;
combined2 = combined;
combined3 = combined;

ax1 = axes;
combined1(:,1:182) = annfp;
combined1(:,183:192) = 1;
combined1(:,193:end) = anntp;
A = imagesc(combined1)
colormap(ax1, cmap)
% colorbar;
A.AlphaData = 0.9;
hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
end

ax2 = axes;
combined2(:,1:182) = valnnfp;
combined2(:,183:192) = 1;
combined2(:,193:end) = valnntp;
B = imagesc(combined2)
colormap(ax2, cmap2)
% colorbar
B.AlphaData = 0.9;

ax3 = axes;
combined3(:,1:182) = MPfp;
combined3(:,183:192) = 1;
combined3(:,193:end) = MPtp;
C = imagesc(combined3)
colormap(ax3, 'gray')
C.AlphaData = 0.55;

linkaxes([ax1,ax2,ax3]);
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];

figure()
imagesc(combined2)
colormap(cmap2)
colorbar
B.AlphaData = 0.6;
%% cronal slice showing difference between MAP18 and MRF-DL
cd T:\Imaging\Multimodal\MRF\Peter
% sfp = 109;
sfp = 111;
cd Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\XXXX\MRF_VBM

cmap = [0 0 0
    0.85 0.225 0.1
    0.85 0.25 0.1
    0.85 0.275 0.1
    0.85 0.3 0.1
    0.85 0.325 0.1
    0.85 0.35 0.1
    0.85 0.375 0.1
    0.9 0.4 0.1
    0.9 0.425 0.1
    0.9 0.45 0.1
    0.9 0.475 0.1
    0.9 0.5 0.1
    0.9 0.525 0.1
    0.9 0.55 0.1
    0.9 0.575 0.1
    0.9 0.6 0.1
    0.9 0.625 0.1
    0.9 0.65 0.1
    0.9 0.675 0.1
    0.9 0.7 0.1]; 

cmap2 = [0 0 0
    1 0.6 0.85
    1 0.675 0.85
    1 0.65 0.85
    1 0.625 0.85
    1 0.6 0.85
    1 0.575 0.85
    1 0.55 0.85
    1 0.525 0.85
    1 0.5 0.85
    1 0.475 0.85
    1 0.45 0.85
    1 0.425 0.85
    1 0.4 0.85
    1 0.375 0.85
    1 0.35 0.85
    1 0.325 0.85
    1 0.3 0.85
    1 0.275 0.85
    1 0.25 0.85
    1 0.225 0.85]; 


MPf = load_untouch_nii('MNI_T1w.nii');
MP = single(MPf.img)*2000;
MPfp = squeeze(MP(:,sfp,:))';
MPfp = flip(MPfp,1);

roif = load_untouch_nii('MNI_ROI_final.nii');
roi = single(roif.img);
roifp = squeeze(roi(:,sfp,:))';
roifp = flip(roifp,1);
% comsize = size(MPfp);
% comsize(2) = comsize(2)*2+10;
% combined = zeros(comsize);
% combined(:,193:end) = roitp;
[B, L] = bwboundaries(roifp,'noholes');

annf = load_untouch_nii('softmax.nii');
ann = single(annf.img);
ann(ann<0.05) = 0;
annfp = squeeze(ann(:,sfp,:))';
annfp = flip(annfp,1);

valnnf = load_untouch_nii('MNI_lesionprob_ANN_fn.nii');
valnn = single(valnnf.img);
valnn(valnn<0.05) = 0;
valnnfp = squeeze(valnn(:,sfp,:))';
valnnfp = flip(valnnfp,1);

figure()
ax1 = axes;
A = imagesc(annfp)
colormap(ax1, cmap)
% colorbar;
A.AlphaData = 0.5;
hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
end

ax2 = axes;
B = imagesc(valnnfp)
colormap(ax2, cmap2)
% colorbar
B.AlphaData = 0.5;


ax3 = axes;
C = imagesc(MPfp)
colormap(ax3, 'gray')
C.AlphaData = 0.55;



linkaxes([ax1,ax2,ax3]);
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];

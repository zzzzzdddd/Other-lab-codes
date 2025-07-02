%% Neural Network Classifier
clear all; close all;
addpath Z:\Imaging\Multimodal\MRF\Peter\MELD_latest_version\matlab
platform_flag = 0; % 0 for windows; 1 for Mac

if platform_flag == 0
    Peter_path = 'Z:\Imaging\Multimodal\MRF\Peter\';
else
    Peter_path = '/Volumes/eegrvw/Imaging/Multimodal/MRF/Peter/';
end

out_dir = fullfile('Z:\Imaging\Multimodal\Myelin', 'gradientsurf');

if ~exist(out_dir)
    mkdir(out_dir)
end
% Directory of patients
subjects_dir = fullfile(Peter_path, 'freesurfer_MELD');

setenv('SUBJECTS_DIR', subjects_dir)
addpath(fullfile(Peter_path, 'MELD_latest_version', 'matlab'));
cd(fullfile(Peter_path, 'MELD_latest_version'));

%change to appropriate prefix
Subset_feature_type = 'SBM'; %select one of {SBM, SBM_FLAIR, SBM_T1overT2, SBM_FLAIR_T1overT2}
% fold = 5;
% time = '1st';
%profix = [num2str(fold) 'fold_' time];
% profix = 'LOO_v5';
profix = 'LOO_v6'; %add the neighbors

subs_all = table2cell(readtable('Subject_gradient.xlsx'))


% load in cortex label left hemi, +1 for freesurfer-matlab indexing
% percentage = 0.5;
Cortex = read_label(['fsaverage_sym'],['lh.cortex']);
Cortex_lesion = Cortex(:,1)+1; % full data 
Cortex_Normal = Cortex_lesion;


All_features={'.Z_by_controls.thickness_z_on_lh.sm10.mgh'; '.Z_by_controls.lh-rh.thickness_z.sm10.mgh';...
    '.Z_by_controls.w-g.pct_z_on_lh.sm10.mgh';'.Z_by_controls.lh-rh.w-g.pct_z.sm10.mgh';...
    '.Z_by_controls.curv_on_lh.mgh';'.Z_by_controls.sulc_on_lh.mgh';...
    '.Z_by_controls.gm_FLAIR_0.75_z_on_lh.sm10.mgh';'.Z_by_controls.gm_FLAIR_0.5_z_on_lh.sm10.mgh';...
    '.Z_by_controls.gm_FLAIR_0.25_z_on_lh.sm10.mgh';'.Z_by_controls.gm_FLAIR_0_z_on_lh.sm10.mgh';...
    '.Z_by_controls.wm_FLAIR_0.5_z_on_lh.sm10.mgh';'.Z_by_controls.wm_FLAIR_1_z_on_lh.sm10.mgh';...
    '.Z_by_controls.lh-rh.gm_FLAIR_0.75_z.sm10.mgh';'.Z_by_controls.lh-rh.gm_FLAIR_0.5_z.sm10.mgh';...
    '.Z_by_controls.lh-rh.gm_FLAIR_0.25_z.sm10.mgh';'.Z_by_controls.lh-rh.gm_FLAIR_0_z.sm10.mgh';...
    '.Z_by_controls.lh-rh.wm_FLAIR_0.5_z.sm10.mgh';'.Z_by_controls.lh-rh.wm_FLAIR_1_z.sm10.mgh';...
    ...
    '.Z_by_controls.gm_MRF_T1_0.75_z_on_lh.sm10.mgh';'.Z_by_controls.gm_MRF_T1_0.5_z_on_lh.sm10.mgh';...
    '.Z_by_controls.gm_MRF_T1_0.25_z_on_lh.sm10.mgh';'.Z_by_controls.gm_MRF_T1_0_z_on_lh.sm10.mgh';...
    '.Z_by_controls.wm_MRF_T1_0.5_z_on_lh.sm10.mgh';'.Z_by_controls.wm_MRF_T1_1_z_on_lh.sm10.mgh';...
    '.Z_by_controls.lh-rh.gm_MRF_T1_0.75_z.sm10.mgh';'.Z_by_controls.lh-rh.gm_MRF_T1_0.5_z.sm10.mgh';...
    '.Z_by_controls.lh-rh.gm_MRF_T1_0.25_z.sm10.mgh';'.Z_by_controls.lh-rh.gm_MRF_T1_0_z.sm10.mgh';...
    '.Z_by_controls.lh-rh.wm_MRF_T1_0.5_z.sm10.mgh';'.Z_by_controls.lh-rh.wm_MRF_T1_1_z.sm10.mgh';...
    ...
    '.Z_by_controls.gm_MRF_T2_0.75_z_on_lh.sm10.mgh';'.Z_by_controls.gm_MRF_T2_0.5_z_on_lh.sm10.mgh';...
    '.Z_by_controls.gm_MRF_T2_0.25_z_on_lh.sm10.mgh';'.Z_by_controls.gm_MRF_T2_0_z_on_lh.sm10.mgh';...
    '.Z_by_controls.wm_MRF_T2_0.5_z_on_lh.sm10.mgh';'.Z_by_controls.wm_MRF_T2_1_z_on_lh.sm10.mgh';...
    '.Z_by_controls.lh-rh.gm_MRF_T2_0.75_z.sm10.mgh';'.Z_by_controls.lh-rh.gm_MRF_T2_0.5_z.sm10.mgh';...
    '.Z_by_controls.lh-rh.gm_MRF_T2_0.25_z.sm10.mgh';'.Z_by_controls.lh-rh.gm_MRF_T2_0_z.sm10.mgh';...
    '.Z_by_controls.lh-rh.wm_MRF_T2_0.5_z.sm10.mgh';'.Z_by_controls.lh-rh.wm_MRF_T2_1_z.sm10.mgh';...
    ...
    '.Z_by_controls.gm_MRF_T2_tr_0.75_z_on_lh.sm10.mgh';'.Z_by_controls.gm_MRF_T2_tr_0.5_z_on_lh.sm10.mgh';...
    '.Z_by_controls.gm_MRF_T2_tr_0.25_z_on_lh.sm10.mgh';'.Z_by_controls.gm_MRF_T2_tr_0_z_on_lh.sm10.mgh';...
    '.Z_by_controls.wm_MRF_T2_tr_0.5_z_on_lh.sm10.mgh';'.Z_by_controls.wm_MRF_T2_tr_1_z_on_lh.sm10.mgh';...
    '.Z_by_controls.lh-rh.gm_MRF_T2_tr_0.75_z.sm10.mgh';'.Z_by_controls.lh-rh.gm_MRF_T2_tr_0.5_z.sm10.mgh';...
    '.Z_by_controls.lh-rh.gm_MRF_T2_tr_0.25_z.sm10.mgh';'.Z_by_controls.lh-rh.gm_MRF_T2_tr_0_z.sm10.mgh';...
    '.Z_by_controls.lh-rh.wm_MRF_T2_tr_0.5_z.sm10.mgh';'.Z_by_controls.lh-rh.wm_MRF_T2_tr_1_z.sm10.mgh';...
    ...
    '.Z_by_controls.lh-rh.pial.K_filtered_2_z.sm20.mgh';'.Z_by_controls.pial.K_filtered_2_z_on_lh.sm20.mgh'
    ...
    '.gm_MRF_T1_0.75_on_lh.sm10.mgh';'.gm_MRF_T1_0.5_on_lh.sm10.mgh';...
    '.gm_MRF_T1_0.25_on_lh.sm10.mgh';'.gm_MRF_T1_0_on_lh.sm10.mgh';...
    '.wm_MRF_T1_0.5_on_lh.sm10.mgh';'.wm_MRF_T1_1_on_lh.sm10.mgh';...
    '.lh-rh.gm_MRF_T1_0.75.sm10.mgh';'.lh-rh.gm_MRF_T1_0.5.sm10.mgh';...
    '.lh-rh.gm_MRF_T1_0.25.sm10.mgh';'.lh-rh.gm_MRF_T1_0.sm10.mgh';...
    '.lh-rh.wm_MRF_T1_0.5.sm10.mgh';'.lh-rh.wm_MRF_T1_1.sm10.mgh';...
    ...
    '.gm_MRF_T2_0.75_on_lh.sm10.mgh';'.gm_MRF_T2_0.5_on_lh.sm10.mgh';...
    '.gm_MRF_T2_0.25_on_lh.sm10.mgh';'.gm_MRF_T2_0_on_lh.sm10.mgh';...
    '.wm_MRF_T2_0.5_on_lh.sm10.mgh';'.wm_MRF_T2_1_on_lh.sm10.mgh';...
    '.lh-rh.gm_MRF_T2_0.75.sm10.mgh';'.lh-rh.gm_MRF_T2_0.5.sm10.mgh';...
    '.lh-rh.gm_MRF_T2_0.25.sm10.mgh';'.lh-rh.gm_MRF_T2_0.sm10.mgh';...
    '.lh-rh.wm_MRF_T2_0.5.sm10.mgh';'.lh-rh.wm_MRF_T2_1.sm10.mgh';...
    ...
    '.gm_MRF_T2_tr_0.75_on_lh.sm10.mgh';'.gm_MRF_T2_tr_0.5_on_lh.sm10.mgh';...
    '.gm_MRF_T2_tr_0.25_on_lh.sm10.mgh';'.gm_MRF_T2_tr_0_on_lh.sm10.mgh';...
    '.wm_MRF_T2_tr_0.5_on_lh.sm10.mgh';'.wm_MRF_T2_tr_1_on_lh.sm10.mgh';...
    '.lh-rh.gm_MRF_T2_tr_0.75.sm10.mgh';'.lh-rh.gm_MRF_T2_tr_0.5.sm10.mgh';...
    '.lh-rh.gm_MRF_T2_tr_0.25.sm10.mgh';'.lh-rh.gm_MRF_T2_tr_0.sm10.mgh';...
    '.lh-rh.wm_MRF_T2_tr_0.5.sm10.mgh';'.lh-rh.wm_MRF_T2_tr_1.sm10.mgh';...
    ...
    '.Z_by_controls.gm_T1_over_T2_0.75_z_on_lh.sm10.mgh';'.Z_by_controls.gm_T1_over_T2_0.5_z_on_lh.sm10.mgh';...
    '.Z_by_controls.gm_T1_over_T2_0.25_z_on_lh.sm10.mgh';'.Z_by_controls.gm_T1_over_T2_0_z_on_lh.sm10.mgh';...
    '.Z_by_controls.wm_T1_over_T2_0.5_z_on_lh.sm10.mgh';'.Z_by_controls.wm_T1_over_T2_1_z_on_lh.sm10.mgh';...
    '.Z_by_controls.lh-rh.gm_T1_over_T2_0.75_z.sm10.mgh';'.Z_by_controls.lh-rh.gm_T1_over_T2_0.5_z.sm10.mgh';...
    '.Z_by_controls.lh-rh.gm_T1_over_T2_0.25_z.sm10.mgh';'.Z_by_controls.lh-rh.gm_T1_over_T2_0_z.sm10.mgh';...
    '.Z_by_controls.lh-rh.wm_T1_over_T2_0.5_z.sm10.mgh';'.Z_by_controls.lh-rh.wm_T1_over_T2_1_z.sm10.mgh';...
    };

%MEASURES
% This is the set of measures you want to include in your classifier
% MRF T2 all use MRF_T2_tr one
SBM = [1:6, 55, 56];
SBM_MRF = [57:68, 81:92];
All = 1:length(All_features);

% lesion_type = [];%original one
lesion_type = '_sm';%smooth one
% lesion_type = '_dil_2';%smooth one
% lesion_type = '_dil_5';%smooth one
% lesion_type = '_dil_10';%smooth one

if strcmp(Subset_feature_type, 'SBM')
    Sets = {'SBM_MRF'}; %SBM_nMRF_v2 is for pyBCE_lr3_SMOTE_dropM and later version
end

cd(subjects_dir)


%% Training Subjects

for SetNumber=1
    Set = eval(Sets{SetNumber});
    SetName = Sets{SetNumber};

    % Normal vertex
    Normal_HC_set = [];
    Measures = All_features(Set);
    NumberOfMeasures = length(Measures);

%     for order = 1:length(HCs_all)
%         NormalSubject = HCs_all{order};
% 
%         if randi(2) == 1
%             h1='lh';
%         else
%             h1='rh';
%         end 
% 
%         %load in all measures for hemisphere
%         for L=1:NumberOfMeasures
%             M = MRIread(['',NormalSubject,'/xhemi/surf/',h1,'',Measures{L},'']);
%             Normal_HC(L,:) = M.vol(Cortex_lesion);
%             clear M 
%         end  
% 
%         Normal_HC_set = [Normal_HC_set, Normal_HC];
%         clear h1 Normal_HC NormalSubject
%     end
%     Normal_HC_set = Normal_HC_set';

    for order = 1:length(subs_all)
        TestSubject = subs_all{order};

        % remove patients
        subset = subs_all;
        Remove = {TestSubject};
        ind = find(ismember(subset, Remove));
        subset(ind) = [];
        clear Remove ind

        disp(['Creating training data for Subject ' TestSubject ' in FT_Set [' SetName ']']);

        if ~exist(fullfile(out_dir, [SetName '_' profix]))
           mkdir(fullfile(out_dir, [SetName '_' profix])) 
        end

        % patient vertex
        for s = 1:length(subset)
            sub = subset(s);
            sub = cell2mat(sub);

            %Set h1 to be lesional hemisphere
            if exist(['',sub,'/xhemi/surf_meld/lh.on_lh.lesion.mgh'])
                h1='lh';
                h2='rh';
            elseif exist(['',sub,'/xhemi/surf_meld/rh.on_lh.lesion.mgh'])
                h1='rh';
                h2='lh';
            else
                display(['error with lesion'])
                break
            end

            %%Load in overlay files from non-lesional hemisphere
            if ~strcmp(sub, 'P83_14516')
                Normal=zeros(length(NumberOfMeasures),length(Cortex_Normal));
                for L=1:NumberOfMeasures
                    M=MRIread(['',sub,'/xhemi/surf/',h2,'',Measures{L},'']);
                    Normal(L,:)=M.vol(Cortex_Normal);
                    clear M
                end
            else
                Normal = [];
            end

            %%Load in lesion label and get indices of vertices
            Lesion = MRIread(['',sub,'/xhemi/surf_meld/',h1,'.on_lh.lesion' lesion_type '.mgh']);
            [a,Lesion,c] = find(Lesion.vol == 1);
            Lesion = intersect(Lesion, Cortex_lesion);
            Lesion = Lesion';

            %%Load in data from lesional hemisphere
            LesionData = zeros(length(NumberOfMeasures),length(Lesion));
            for L = 1:NumberOfMeasures
                M = MRIread(['',sub,'/xhemi/surf/',h1,'',Measures{L},'']);
                LesionData(L,:) = M.vol(Lesion);
                visual = horzcat(, Cortex(:, 2:4), LesionData')
                visual = pointCloud(Cortex(:, 2:4));
                pcshow(visual)
                clear M
            end

            %%Load all data together as a single matrix
%             Combined = [Normal, LesionData];
%             Combined = Combined';

            %%Add Normal and Lesional data to data of all other subjects
%             Multi = [Multi; Combined];

%             Binary(1:size(Normal, 2)) = 1;
%             Binary(size(Normal, 2)+1:size(Normal, 2)+size(LesionData, 2)) = 0;
%             Binary = Binary';
%             Score = [Score; Binary(:)];

            clear Binary Combined Normal Lesion LesionData rand_idx sub
        end

        % Use principal component analysis to pick number of nodes in neural
        % network - number of components explaining 99% of the variance
%             [COEFF,SCORE,latent,tsquare] = pca(Multi);
%             PCA=cumsum(latent)./sum(latent);  
%             hiddenLayerSize=[length(find(PCA < 0.995)), length(find(PCA < 0.92)), length(find(PCA < 0.68))];
%         rand_idx = randi(size(Normal_HC_set, 1), size(Multi, 1), 1);
%         Multi = cat(1, Multi, Normal_HC_set(rand_idx, :));
%         Score = cat(1, Score, ones(size(rand_idx, 1), 1));
%         Score(:,2) = Score(:,1).*-1+1;

%         save(fullfile(out_dir, [SetName '_' profix], [TestSubject lesion_type '_data_Peter.mat']), 'X_train', 'Y_train', '-v7.3')
%         clear X_train Y_train TestSubject subset
    end
%     clear Normal_HC_set Measures NumberOfMeasures SetName Set
end

%% HC Subjects
for order = 1:length(HCs_all)
    TestSubject = HCs_all{order};
    
    for SetNumber = 1:length(Sets)
        Set = eval(Sets{SetNumber});
        SetName = Sets{SetNumber};
        disp(['Creating HC data for Subject ' TestSubject ' in FT_Set [' SetName ']']);
        if ~exist(fullfile(out_dir, [SetName '_' profix '_HC']))
           mkdir(fullfile(out_dir, [SetName '_' profix '_HC'])) 
        end
        %%selecting measures
        Measures = All_features(Set);
        NumberOfMeasures = length(Measures);
        
        for hemi=1:2
            %run on both hemispheres
            if hemi==1
                h1='lh';
            else
                h1='rh';
            end 

            %load in all measures for hemisphere
            for L = 1:NumberOfMeasures
                M = MRIread(['',TestSubject,'/xhemi/surf/',h1,'',Measures{L},'']);
                Normal(L,:) = M.vol(Cortex_lesion);
            end  
            eval(['X_test_' h1 ' = Normal'';']);    
            
            clear X Normal h1
        end
        Y_test_lh(:, 2) = zeros([size(X_test_lh, 1), 1]);
        Y_test_lh(:, 1) = ones([size(X_test_lh, 1), 1]);
        Y_test_rh(:, 2) = zeros([size(X_test_rh, 1), 1]);
        Y_test_rh(:, 1) = ones([size(X_test_rh, 1), 1]);
                
        save(fullfile(out_dir, [SetName '_' profix '_HC'], [TestSubject lesion_type '_data.mat']), 'X_test_lh', 'X_test_rh', 'Y_test_lh', 'Y_test_rh', '-v7.3');
        clear Y_test_lh Y_test_rh Lesion_lh Lesion_rh
        clear Measures NumberOfMeasures SetName Set
    end
    clear TestSubject
end
%% x-fold cross validation
for Fcn_order = 1
    for order = 1:length(Subject_folds)
        TestSubjects = Subject_folds{order};
        
        % remove patients
        subset = subs_all;
        Remove = TestSubjects;
        ind = find(ismember(subset,Remove));
        subset(ind)=[];
        clear Remove ind
        
        for SetNumber=1:length(Sets)
            Set = eval(Sets{SetNumber});
            SetName = Sets{SetNumber};
            disp(['Fold_' num2str(order)  '   ' TrainFcn_set{Fcn_order} '   ' SetName]);
            
            if ~exist(fullfile(out_dir, [SetName '_' num2str(fold) 'fold_' time]))
                mkdir(fullfile(out_dir, [SetName '_' num2str(fold) 'fold_' time])) 
            end
            
            %%selecting measures
            Measures=All_features(Set);
            NumberOfMeasures=length(Measures);
            
            Multi=[];
            Score=[];
            
            % training data
            for s = 1:length(subset)
                sub = subset(s);
                sub = cell2mat(sub);

                %Set h1 to be lesional hemisphere
                if exist(['',sub,'/xhemi/surf_meld/lh.on_lh.lesion.mgh'])
                    h1='lh';
                    h2='rh';
                elseif exist(['',sub,'/xhemi/surf_meld/rh.on_lh.lesion.mgh'])
                    h1='rh';
                    h2='lh';
                else
                    display(['error with lesion'])
                    break
                end

                %%Load in overlay files from non-lesional hemisphere
                if ~strcmp(sub, 'P83_14516')
                    Normal=zeros(length(NumberOfMeasures),length(Cortex_Normal));
                    for L=1:NumberOfMeasures
                        M=MRIread(['',sub,'/xhemi/surf/',h2,'',Measures{L},'']);
                        Normal(L,:)=M.vol(Cortex_Normal);
                        clear M
                    end
                else
                    Normal = [];
                end

                %%Load in lesion label and get indices of vertices
                Lesion=MRIread(['',sub,'/xhemi/surf_meld/',h1,'.on_lh.lesion' lesion_type '.mgh']);
                [a,Lesion,c]=find(Lesion.vol==1);
                Lesion=intersect(Lesion,Cortex_lesion);
                Lesion=Lesion';

                %%Load in data from lesional hemisphere
                LesionData=zeros(length(NumberOfMeasures),length(Lesion));
                for L=1:NumberOfMeasures
                    M=MRIread(['',sub,'/xhemi/surf/',h1,'',Measures{L},'']);
                    LesionData(L,:)=M.vol(Lesion);
                    clear M
                end

                %%Load all data together as a single matrix
                Combined=[Normal,LesionData];
                Combined=Combined';

                %%Add Normal and Lesional data to data of all other subjects
                Multi=[Multi;Combined];

                Binary(1:size(Normal, 2))=1;
                Binary(size(Normal, 2)+1:size(Normal, 2)+size(LesionData, 2))=0;
                Binary=Binary';
                Score=[Score; Binary(:)];

                clear Binary Combined Normal Lesion LesionData rand_idx sub
            end
            
            Score(:,2)=Score(:,1).*-1+1; 
            
            X_train = Multi;
            Y_train = Score;
            save(fullfile(out_dir, [SetName '_' num2str(fold) 'fold_' time], ['Fold' num2str(order) lesion_type '_data.mat']), 'X_train', 'Y_train', '-v7.3')
            clear X_train Y_train
            
            %testing data
            for s = 1:length(TestSubjects)
                sub = TestSubjects{s};
                
                %Subjects list minus test subject
                for hemi=1:2
                    %run on both hemispheres
                    if hemi==1
                        h1='lh';
                        h2='rh';
                    else
                        h1='rh';
                        h2='lh';
                    end 

                    %load in all measures for hemisphere
                    for L = 1:NumberOfMeasures
                        M = MRIread(['', sub,'/xhemi/surf/',h1,'',Measures{L},'']);
                        Normal(L,:)=M.vol(Cortex_lesion);
                    end

                    eval(['X_test_' h1 ' = Normal'';']);
                    X{1}=Normal;

                    clear M Y Out X Normal
                end
                
                if exist(['', sub,'/xhemi/surf_meld/lh.on_lh.lesion.mgh'])
                    Lesion_lh=MRIread(['',sub,'/xhemi/surf_meld/','lh','.on_lh.lesion' lesion_type '.mgh']);
                    Y_test_lh(:, 2) = Lesion_lh.vol(Cortex_lesion)';
                    Y_test_lh(:, 1) = Y_test_lh(:, 2)*-1+1;
                    Y_test_rh(:, 1) = ones(size(Y_test_lh(:, 1)));
                    Y_test_rh(:, 2) = zeros(size(Y_test_lh(:, 2)));
                else
                    Lesion_rh=MRIread(['',sub,'/xhemi/surf_meld/','rh','.on_lh.lesion' lesion_type '.mgh']);
                    Y_test_rh(:, 2) = Lesion_rh.vol(Cortex_lesion)';
                    Y_test_rh(:, 1) = Y_test_rh(:, 2)*-1+1;
                    Y_test_lh(:, 1) = ones(size(Y_test_rh(:, 1)));
                    Y_test_lh(:, 2) = zeros(size(Y_test_rh(:, 2)));
                end
                
                save(fullfile(out_dir, [SetName '_' num2str(fold) 'fold_' time], ['Fold' num2str(order) '_TestSub_' sub lesion_type '_data.mat']), 'X_test_lh', 'X_test_rh', 'Y_test_lh', 'Y_test_rh', '-v7.3')

                clear X_test_lh X_test_rh Y_test_lh Y_test_rh Lesion_lh Lesion_rh
            end
            clear Set SetName Measures NumberOfMeasures
        end
        clear subset TestSubjects
    end
end


%% HCs training data generation 
for order = 1:length(HCs_all)
    TestSubject = HCs_all{order};

    for SetNumber=1:length(Sets)
        Set=eval(Sets{SetNumber});
        SetName=Sets{SetNumber};
        disp([TestSubject '   ' SetName]);

        if ~exist(fullfile(out_dir, SetName))
            mkdir(fullfile(out_dir, SetName));
        end

        %%selecting measures
        Measures = All_features(Set);
        NumberOfMeasures = length(Measures);
        Header_type = strip(erase(Measures, ["Z_by_controls", "mgh"]), 'both', '.');
        
        %Subjects list minus test subject
        for hemi=1:2
            %run on both hemispheres
            if hemi==1
                h1='lh';
            else
                h1='rh';
            end 

            %load in all measures for hemisphere
            for L=1:NumberOfMeasures
                file_H1 = dir(fullfile(TestSubject, 'xhemi', 'surf', [h1 '.*' Header_type{L} '*.mgh']));
                M = MRIread(fullfile(TestSubject, 'xhemi', 'surf', file_H1.name));
        
                Normal(L,:) = M.vol(Cortex_lesion);
                clear file_H1 M
            end

            eval(['X_test_' h1 ' = Normal'';']);
            clear Normal
        end
        Y_test_lh(:, 2) = zeros([size(X_test_lh, 1), 1]);
        Y_test_lh(:, 1) = ones([size(X_test_lh, 1), 1]);
        Y_test_rh(:, 2) = zeros([size(X_test_rh, 1), 1]);
        Y_test_rh(:, 1) = ones([size(X_test_rh, 1), 1]);
        
        save(fullfile(out_dir, SetName, [TestSubject lesion_type '_data.mat']), 'X_test_lh', 'X_test_rh', 'Y_test_lh', 'Y_test_rh', '-v7.3');
        
        clear SetName Set Measures NumberOfMeasures
        clear X_test_lh X_test_rh Y_test_lh Y_test_rh Header_type
    end
    clear TestSubject
end

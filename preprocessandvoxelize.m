function preprocessandvoxelize(destination,source)
%% Load data
volLoc = [source filesep 'imagesTr'];
lblLoc = [source filesep 'labelsTr'];

% If the directory for preprocessed data does not exist, or only a partial
% set of the data files have been processed, process the data.
if ~exist(destination,'dir') || proceedWithPreprocessing(destination)

    mkdir(fullfile(destination,'imagesTr'));
    mkdir(fullfile(destination,'labelsTr'));
    
    mkdir(fullfile(destination,'imagesVal'));
    mkdir(fullfile(destination,'labelsVal'));
    
    mkdir(fullfile(destination,'imagesTest'));
    mkdir(fullfile(destination,'labelsTest'));
    
    labelReader = @(x) (niftiread(x) > 0);
    volReader = @(x) niftiread(x);
    volds = imageDatastore(volLoc, ...
        'FileExtensions','.nii','ReadFcn',volReader);
    classNames = ["normal","lesion"];
    pixelLabelID = [0 1];
    pxds = pixelLabelDatastore(lblLoc,classNames, pixelLabelID, ...
        'FileExtensions','.nii','ReadFcn',labelReader);
    reset(volds);
    reset(pxds);
    
    %% Crop relevant region
    NumFiles = length(pxds.Files);
    id = 1;
    while hasdata(pxds)
        outL = readNumeric(pxds);
        outV = read(volds);

        incropVol = outV(17:164,17:202,1:157,:);
        mask = incropVol == 0;
        cropVol = channelWisePreProcess(incropVol);
        
        % Set the nonbrain region to 0.
        cropVol(mask) = 0;
        cropLabel = outL(17:164,17:202,1:157,:);
        
        % Split data into training, validation and test sets
%         if (id < floor(0.3*NumFiles))|(id > floor(0.5*NumFiles))
%             imDir    = fullfile(destination,'imagesTr','Sub');
%             labelDir = fullfile(destination,'labelsTr','Sub');
%         else 
%             imDir    = fullfile(destination,'imagesTest','Sub');
%             labelDir = fullfile(destination,'labelsTest','Sub');


        imDir    = fullfile(destination,'imagesTr','Sub');
        labelDir = fullfile(destination,'labelsTr','Sub'); 
        save([imDir num2str(id,'%.3d') '.mat'],'cropVol');
        save([labelDir num2str(id,'%.3d') '.mat'],'cropLabel');
        imDir    = fullfile(destination,'imagesTest','Sub');
        labelDir = fullfile(destination,'labelsTest','Sub');
        save([imDir num2str(id,'%.3d') '.mat'],'cropVol');
        save([labelDir num2str(id,'%.3d') '.mat'],'cropLabel');

%         elseif (id < floor(0.25*NumFiles))|(id < floor(0.55*NumFiles))
%             imDir    = fullfile(destination,'imagesTest','Sub');
%             labelDir = fullfile(destination,'labelsTest','Sub');
%         else
%             imDir    = fullfile(destination,'imagesVal','Sub');
%             labelDir = fullfile(destination,'labelsVal','Sub');
%         end
        

%         save([imDir num2str(id,'%.3d') '.mat'],'cropVol');
%         save([labelDir num2str(id,'%.3d') '.mat'],'cropLabel');


        id=id+1;       
   end
end
% cnt
end

function out = channelWisePreProcess(in)
% 
out = in;
% chn_Mean = mean(in(in~=0),[1 2 3]);
% chn_Std = std(in(in~=0),0,[1 2 3]);
% out = (in - chn_Mean)./chn_Std;

rangeMin = -8;
rangeMax = 8;
morphs = out(:,:,:,3:8);
morphs(morphs > rangeMax) = rangeMax;
morphs(morphs < rangeMin) = rangeMin;
out(:,:,:,3:8) = morphs;


% out = (out - rangeMin) / (rangeMax - rangeMin);

% Rescale the data to the range [0, 1].
minchannels = min(out,[],[1 2 3]);
maxchannels = max(out,[],[1 2 3]);
minchannels(:,:,:,3:8) = -8;
maxchannels(:,:,:,3:8) = 8;
out = (out - minchannels) ./ (maxchannels - minchannels);
end

function out = proceedWithPreprocessing(destination)
totalNumFiles = 162;
numFiles = 0;
if exist(fullfile(destination,'imagesTr'),'dir')
    tmp1 = dir(fullfile(destination,'imagesTr'));
    numFiles = numFiles + sum(~vertcat(tmp1.isdir));
end

if exist(fullfile(destination,'imagesVal'),'dir')
    tmp1 = dir(fullfile(destination,'imagesVal'));
    numFiles = numFiles + sum(~vertcat(tmp1.isdir));
end

if exist(fullfile(destination,'imagesTest'),'dir')
    tmp1 = dir(fullfile(destination,'imagesTest'));
    numFiles = numFiles + sum(~vertcat(tmp1.isdir));
end

% If total number of preprocessed files is not equal to the number of
% files in the dataset, perform preprocessing. Otherwise, preprocessing has
% already been completed and can be skipped.
out = (numFiles ~= totalNumFiles);
end
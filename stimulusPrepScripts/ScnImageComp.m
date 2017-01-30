% Script to find most spatially dissimilar scene images on the categories
% dataset.
clearvars;
clc;

masterDir   = '~/Google Drive/Research/tACS/tACS_ER_task/stim/';
origStimDir   = '~/Google Drive/Research/tACS/tACS_ER_task/stim/scene_categories_orig/';
newStimDir   = '~/Google Drive/Research/tACS/tACS_ER_task/stim/scene_categories/';

SceneCat    = [];
% Indoor
SceneCat{1} = {'office','kitchen','store','bedroom','livingroom'};
% Outdoor
SceneCat{2} = {'suburb','highway','opencountry','industrial','coast','insidecity',...
    'street','forest','mountain','tallbuilding'};
nSceneCat   = cellfun(@numel,SceneCat);
SceneCatNames = {'indoor','outdoor'};
%
ImgNames        = cell(2,1);
NumImages       = cell(2,1);
StimDirsOrig    = cell(2,1);
StimDirsNew     = cell(2,1);
ImgSizes        = cell(2,1);
ImgRescale      = 256;
for ii = 1:2
    for jj = 1:nSceneCat(ii)
        StimDirsOrig{ii,jj} = strcat(origStimDir,SceneCatNames{ii},'/',SceneCat{ii}{jj},'/');
        StimDirsNew{ii,jj}  = strcat(newStimDir,SceneCatNames{ii},'/',SceneCat{ii}{jj},'/');
        
        temp = dir([StimDirsOrig{ii,jj} 'image_*.jpg']);
        ImgNames{ii}{jj} = {temp.name}';
        NumImages{ii}(jj) = numel(temp);
        nImg = NumImages{ii}(jj);
        for kk = 1:nImg
            temp = imfinfo([StimDirsOrig{ii,jj} ImgNames{ii}{jj}{kk}]);
            ImgSizes{ii}{jj}(kk,:) = [temp.Width temp.Height];
        end
    end
end

%% crop to minimum size in that category, re-scale to 256, and create new files.
if 0
    for ii = 1:2
        for jj = 1:nSceneCat(ii)
            nImg = NumImages{ii}(jj);
            minImgSize = min(ImgSizes{ii}{jj});
            t=tic;
            for kk = 1:nImg
                img1S   = ImgSizes{ii}{jj}(kk,:);
                pad     = round((img1S-minImgSize)/2);
                % only crop/resize for images that exceed the size requirements
                if (sum(pad)>0) || any(ImgRescale ~= minImgSize)
                    img1     = imread([StimDirsOrig{ii,jj} ImgNames{ii}{jj}{kk}]);
                    img2     = imcrop(img1, [pad,minImgSize]);
                    img3     = imresize(img2, [ImgRescale ImgRescale],'method','lanczos3');
                    img3     = uint8(img3(:,:,1));
                    imwrite(img3,[StimDirsNew{ii,jj} ImgNames{ii}{jj}{kk}]);
                end
            end
            fprintf('\ntime to resize and crop category %s = %g secs \n',SceneCat{ii}{jj},toc(t))
        end
    end
end

%% manual delete of images
% 1) images w people or animals in them
% 2) images w readable signs
% 3) images w poor quality
% 4) images w repeats or take from a different angle

%% Compute correlations to find spatial similarity within stimulus category
% to do on cropped images
ImgCorr         = cell(2,1);
for ii = 1:2
    for jj = 1:nSceneCat(ii)
        temp = dir([StimDirsNew{ii,jj} 'image_*.jpg']);
        ImgNames{ii}{jj} = {temp.name}';
        
        NumImages{ii}(jj) = numel(temp);
        nImg = NumImages{ii}(jj);
        
        C  = nan(nImg);
        t=tic;       
        for kk1 = 1:(nImg-1)
            img1 = imread([StimDirsNew{ii,jj} ImgNames{ii}{jj}{kk1}]);
            for kk2 = kk1+1:nImg
                img2 = imread([StimDirsNew{ii,jj} ImgNames{ii}{jj}{kk2}]);
                C(kk1,kk2) = corr2(img1,img2);
            end
        end
        fprintf('\ntime to create category %s = %g secs \n',SceneCat{ii}{jj},toc(t))
        ImgCorr{ii}{jj} = C;
    end
end
save([newStimDir,'imageCorrsv2'],'ImgCorr','ImgNames','SceneCat')
%% Get image sampling weights based on correlation and
corr_thr = 0.4;
W = cell(2,1); % weights
R = cell(2,1); % ranks
for ii = 1:2
    for jj = 1:nSceneCat(ii)
        C = ImgCorr{ii}{jj};  % get correlation scores
        C(isnan(C))=0;
        binC = abs((C+C'))>corr_thr; % obtain the number of images that each
        sumBC = sum(binC);
        xx=1./(1+sumBC);
        W{ii}{jj} = xx/sum(xx); % get probability sampling weights
        [~,R{ii}{jj}] = sort(sumBC);
    end
end
save([newStimDir,'imageSamplingWeightsv2'],'W','ImgNames','SceneCat')
save([newStimDir,'imageRanksv2'],'R','ImgNames','SceneCat')

%% Create new dir for in vs outdoor cat. % behav_v13 version
clear all;
masterDir   = '~/Google Drive/Research/tACS/tACS_ER_task/stim/';
load([masterDir 'scene_categories/imageRanksv2.mat'])
allStimDir      = fullfile(masterDir,'scene_categories');

nStimPerCat     = 300;
nStimPerSubCat  = 60;
selSceneCats = [1:5;1 2 3 5 9]; % subselection of categories
SceneCatNames = {'indoor','outdoor'};

ImgMat = zeros(2,nStimPerCat,256,256,'uint8');
selImgNames = cell(2,nStimPerCat);

for ii = 1:2
    cnt = 1;
    for jj2 = 1:5
        jj = selSceneCats(ii,jj2);
        catDir = fullfile(allStimDir,SceneCatNames{ii},SceneCat{ii}{jj});
        for kk =1:nStimPerSubCat
            imgID = R{ii}{jj}(kk);
            imgName = ImgNames{ii}{jj}{imgID};
            ImgMat(ii,cnt,:,:) = imread([catDir '/' imgName]);
            selImgNames{ii,cnt} = strcat(SceneCat{ii}{jj},'_',imgName);
            cnt = cnt+1;
        end
    end
end
save([allStimDir,'/selInOutImgsV2'],'ImgMat','selImgNames');

%% Create new file for images for behav_v14 version. man made vs natural
clear all;
masterDir   = '~/Google Drive/Research/tACS/tACS_ER_task/stim/';
load([masterDir 'scene_categories/imageRanksv2.mat'])
allStimDir      = fullfile(masterDir,'scene_categories');

nStimPerCat     = 300;
nStimPerSubCat  = 100;
selSceneCats = [1 2 4; 5 8 9]; % subselection of categories
SceneCatNames = {'outdoor','outdoor'}; % man made vs natural
%SceneCatNames = {'mm','na'}; % man made vs natural

ImgMat = zeros(2,nStimPerCat,256,256,'uint8');
selImgNames = cell(2,nStimPerCat);

for ii = 1:2
    cnt = 1;
    for jj2 = 1:3
        jj = selSceneCats(ii,jj2);
        catDir = fullfile(allStimDir,'outdoor',SceneCat{2}{jj});
        for kk =1:nStimPerSubCat
            imgID = R{2}{jj}(kk);
            imgName = ImgNames{2}{jj}{imgID};
            ImgMat(ii,cnt,:,:) = imread([catDir '/' imgName]);
            selImgNames{ii,cnt} = strcat(SceneCat{2}{jj},'_',imgName);
            cnt = cnt+1;
        end
    end
end
save([allStimDir,'/selInOutImgs_v14'],'ImgMat','selImgNames');
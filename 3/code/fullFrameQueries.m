% Computer Vision Project 3
% Turned in by <Michael Single>
% Legi: 08-917-445


% flush context
clc
clear all;
close all;

% load dependencies
run lib/vlfeat-0.9.19/toolbox/vl_setup 

addpath('util/');
addpath('lib/');
addpath('src/');

% should we use the breaking bad data set?
useBB = false;

% assumes that such a file exists, can be generated by running
% the function #computeSiftDataOf(..). See rawDescriptorMatches for its
% usage.
wordsPerCluster = 1500;

frameBaseName = 'video_test';
if useBB
    frameBaseName = 'breakingbad2';
end
allSiftFeaturesBase = strcat('data/',frameBaseName,'_all.mat');

load(allSiftFeaturesBase, 'allFrameNames', 'descrCount', 'allPos', ...
    'allOrients', 'allScales', 'allDescr');
%%
[membership, means, rms] = kmeansML(wordsPerCluster, allDescr');
disp('Clustering data has been computed.');

%%
frameCount = length(descrCount);
histData = computeFrameHisto(descrCount, membership, frameCount, wordsPerCluster);
normalizationFScale = sqrt(sum(histData.^2, 2));
randFrames = randperm(frameCount, 3);
for frameIdx = randFrames
    histOfRanFrame = histData(frameIdx, :);
    
    % similarity is dot(hist_sel, hist_other) / (norm(hist_sel)*norm(hist_other))
    normalFrameF = (repmat(norm(histOfRanFrame),frameCount,1).*normalizationFScale);
    scores = dot(repmat(histOfRanFrame, frameCount,1), histData, 2);
    
    % order scores descending, i.e. lower indices refer to better scores.
    [scores, idxs] = sort(scores ./ normalFrameF, 'descend');
    figure('name', 'Full Frame Queries')
    
    currSiftMatFile = strcat(frameBaseName, '_', num2str(idxs(1)), '.mat');
    [currentFrameName, featureCount, positions, ...
        orients, scales, descriptors] = loadOwnSiftDataOf(currSiftMatFile);
    img = imread(strcat('frames/',currentFrameName));
    subplot(2,3,1);
    imshow(img);
    
    tillIdx = regexpi(currentFrameName, '_');
    fraName = currentFrameName(tillIdx+1:end);
    tillIdx = regexpi(fraName, '.png');
    fraName = fraName(1:tillIdx-1);
    title(strcat('Selected Image ', fraName));
    for k=2:6
        currSiftMatFile = strcat(frameBaseName, '_', num2str(idxs(k)), '.mat');
        [currentFrameName, featureCount, positions, ...
            orients, scales, descriptors] = loadOwnSiftDataOf(currSiftMatFile);
        img = imread(strcat('frames/',currentFrameName));
        subplot(2,3,k);
        imshow(img);
        
        tillIdx = regexpi(currentFrameName, '_');
        fraName = currentFrameName(tillIdx+1:end);
        tillIdx = regexpi(fraName, '.png');
        fraName = fraName(1:tillIdx-1);
        title(strcat('Matching Imgage ', fraName, ' with score: ',num2str(scores(k))));
    end

end

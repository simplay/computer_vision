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


useBB = false;

% assumes that such a file exists, can be generated by running
% the function #computeSiftDataOf(..). See rawDescriptorMatches for its
% usage.
allSiftFeaturesBase = 'data/video_test_all.mat';
if useBB
    allSiftFeaturesBase = 'data/breakingbad2_all.mat';
end

wordsPerCluster = 1500;

load(allSiftFeaturesBase, 'allFrameNames', 'descrCount', 'allPos', ...
    'allOrients', 'allScales', 'allDescr');

%% compute clustering data
[membership, means, rms] = kmeansML(wordsPerCluster, allDescr');

%% find vocabulary
randWordIndices = randperm(wordsPerCluster,2);
for wordIdx = randWordIndices
    figure('name', strcat('showing voca of word', num2str(wordIdx)))
    foundFeatureIdxs = find(membership==wordIdx);
    showPatchCount = min(length(foundFeatureIdxs), 25);
    for k=1:showPatchCount
        featureIdx = foundFeatureIdxs(k);
        imgPatchName = strcat('frames/',allFrameNames{featureIdx,:});
        img = rgb2gray(imread(imgPatchName));
        patch = getPatchFromSIFTParameters(allPos(featureIdx,:), allScales(featureIdx), allOrients(featureIdx), img);
        subplot(5,5,k)
        imshow(patch);
    end
end

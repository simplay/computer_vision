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

%% specify initial runtime paramters

% specify data sets
useTwoFrameData = false;
useBreakingBad = true;
shouldComputeDataSets = false;

% init file loading setting.
fileIdx1 = 1;
fileIdx2 = 40; % for e.g. video_test
videoFileName = 'video_test.mp4';
if useBreakingBad
    videoFileName = 'breakingbad2.mp4';
    fileIdx1 = 1;
    fileIdx2 = 140; % select the white tubus in the back
end
tillIdx = regexpi(videoFileName, '\.');
frameBaseName = videoFileName(1:tillIdx-1);

%% compute data sets
if shouldComputeDataSets
    % extracts frames from a given videos and stores them as png files in the
    % folder 'frames'.
    extractVideoFrames(videoFileName)

    % compute sift data for each frame associated to the given video file.
    % stores it as a mat file called 'frameBaseName'.
    computeSiftDataOf(frameBaseName);
end


%% load correspinding sift data and left and right img frames according to specified data-set
if (useTwoFrameData)
    load('twoFrameData.mat');
else
    leftSiftMatFile = strcat(frameBaseName, '_', num2str(fileIdx1), '.mat');
    rightSiftMatFile = strcat(frameBaseName, '_', num2str(fileIdx2), '.mat');
    
    % left frame sift data
    [currentFrameName1, featureCount1, positions1, ...
        orients1, scales1, descriptors1] = loadOwnSiftDataOf(leftSiftMatFile);
    im1 = imread(strcat('frames/',currentFrameName1));
    
    % right frame sift data
    [currentFrameName2, featureCount2, positions2, ...
        orients2, scales2, descriptors2] = loadOwnSiftDataOf(rightSiftMatFile);
    im2 = imread(strcat('frames/',currentFrameName2));
end

% let user specify a region of interest
maskedIndices = selectRegion(im1, positions1);
maskedIdxDescr1 = descriptors1(maskedIndices, :);

% find the indices between significant, closely matching points.
matches = getMatchingFeaturePoints(maskedIdxDescr1, descriptors2);

%% show matches
figure('name', 'showing matching features');
subplot(1,2,1);
imshow(im1);
[pos, orients, scales] = getSelMatchings(matches(:,1), positions1, ...
                                         orients1, scales1, maskedIndices);
% display selection features
displaySIFTPatches(pos, scales, orients, im1);

% display matching features
subplot(1,2,2);
imshow(im2);

% find good matches
[pos, orients, scales] = getSelMatchings(matches(:,2), positions2, orients2, scales2);
displaySIFTPatches(pos, scales, orients, im2);

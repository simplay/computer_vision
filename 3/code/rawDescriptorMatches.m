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

% specify initial runtime paramters
usetTwoFrameData = false;
videoFileName = 'video_test.mp4';
fileIdx1 = 1;
fileIdx2 = 244;

% extracts frames from a given videos and stores them as png files in the
% folder 'frames'.
extractVideoFrames(videoFileName)

% compute sift data for each frame associated to the given video file.
% stores it as a mat file called 'frameBaseName'.
tillIdx = regexpi(videoFileName, '\.');
frameBaseName = videoFileName(1:tillIdx-1);
computeSiftDataOf(frameBaseName);

%%
if (usetTwoFrameData)
    load('twoFrameData.mat');
else
    leftSiftMatFile = strcat(frameBaseName, '_', num2str(fileIdx1), '.mat');
    rightSiftMatFile = strcat(frameBaseName, '_', num2str(fileIdx2), '.mat');

    [currentFrameName1, featureCount1, positions1, ...
        orients1, scales1, descriptors1] = loadOwnSiftDataOf(leftSiftMatFile);
    im1 = imread(strcat('frames/',currentFrameName1));
    
    [currentFrameName2, featureCount2, positions2, ...
        orients2, scales2, descriptors2] = loadOwnSiftDataOf(rightSiftMatFile);
    im2 = imread(strcat('frames/',currentFrameName2));
end


maskedIndices = selectRegion(im1, positions1);
maskedIdxDescr1 = descriptors1(maskedIndices, :);
matches = getMatchingFeaturePoints(maskedIdxDescr1, descriptors2);

%% show matches
figure('name', 'showing matching features');
subplot(1,2,1);
imshow(im1);
[pos, orients, scales] = getSelMatchings(matches(:,1), positions1, ...
                                         orients1, scales1, maskedIndices);
displaySIFTPatches(pos, scales, orients, im1);

subplot(1,2,2);
imshow(im2);
[pos, orients, scales] = getSelMatchings(matches(:,2), positions2, orients2, scales2);
displaySIFTPatches(pos, scales, orients, im2);

clear all
%close all
clc

maxSize = 500;
% We assume images have the same size
left = mean(single(imread('pics/left.jpg')),3);
[M, N] = size(left);
ratio = max(M,N)/maxSize;
left = imresize(left,1/ratio);

right = mean(single(imread('pics/right.jpg')),3);
right = imresize(right,1/ratio);

[rows, cols] = size(left);

figure(1)
clf;
subplot(1,2,1);
imagesc(left);
colormap(gray) 
title('Left Image');

subplot(1,2,2);
imagesc(right);
colormap(gray)
title('Right Image');

% Number of corresponding points.
numPoints = 8;
% Set to true to save the points in the file savedPoints.mat.
savePoints = true;
% Set to true to select new points even if a the file savedPoints.mat
% exists.
selectNewPoints  = false;

% Corresponding points in homogeneous coordinates.
leftPoints = zeros(3, numPoints);
rightPoints = zeros(3, numPoints);

if exist('savedPoints.mat','file') && ~selectNewPoints
    load('savedPoints.mat')
else
    for i=1:numPoints
        disp(['Selection of Position ', num2str(i), ' of ', num2str(numPoints), ':']);
        disp('Mark any Position in the Left image.');
        [x, y] = ginput(1);
        leftPoints(:,i) = [x, y, 1]';
        
        disp('Mark the corresponding Position in the Right image.');
        [x, y] = ginput(1);
        rightPoints(:,i) = [x, y, 1]';
        clc
    end
    
    if savePoints
        save('savedPoints.mat','leftPoints','rightPoints');
    end
end


% Computing the fundamental matrix
F = eightPointsAlgorithm(leftPoints,rightPoints);  



disp('Select a point in the left image. Press ESC to exit.');
while true
    

    [x, y, key] = ginput(1);
    if key==27
        break;
    end
    
    % TODO: Question 4
end

% Compute epipoles
% TODO: Question 5

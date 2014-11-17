clear all
close all
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

figure('Position',[100,100,1024,800])
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
        
        % 2D homogenous coordinates: Positions have a 1 as their last component.
        leftPoints(:,i) = [x, y, 1]';
        
        disp('Mark the corresponding Position in the Right image.');
        [x, y] = ginput(1);
        
        % 2D homogenous coordinates: Positions have a 1 as their last component. 
        rightPoints(:,i) = [x, y, 1]';
        clc
    end
    
    if savePoints
        save('savedPoints.mat','leftPoints','rightPoints');
    end
end


% Computing the fundamental matrix
F = eightPointsAlgorithm(leftPoints,rightPoints); 

% epipolar line left image is equal to the kernel of F, since F*e = 0
% epipolar line right image is equal to the kernel of F', since F'*e' = 0
epipoleLeft = null(F);
epipoleRight = null(F');

% In order to retrieve image coordiantes of epipole positions 
% perform a homogeneous division (e.x/e.w, e.y/e.w).
leftEpipolePosition = epipoleLeft/epipoleLeft(3);
rightEpipolePosition = epipoleRight/epipoleRight(3);

disp('Select a point in the left image. Press ESC to exit.');
while true
    
    % coordinates of points selected in the left image
    [px, py, key] = ginput(1);
    leftPoint = [px; py; 1];
    if key==27
        break;
    end
    
    subplot(1,2,1);
    imagesc(left);
    hold on
    
    % plot clicked point
    plot(px, py,'c.');
    
    % general Form of a line ax + by + c = 0
    % thus y = (-c - ax)/b
    leftEpipolarLineCoeff = cross(epipoleLeft, leftPoint);
    a = leftEpipolarLineCoeff(1); 
    b = leftEpipolarLineCoeff(2); 
    c = leftEpipolarLineCoeff(3);
    
    % left epiline 
    xLeft = [0; size(left,2)]; yLeft = (-c-a*xLeft)/b;
    
    % plot epipolar line in the left image
    plot(xLeft, yLeft, 'r-');
    
    subplot(1,2,2);
    imagesc(right);
    hold on
    
    % general Form of a line ax + by + c = 0
    % thus y = (-c - ax)/b
    rightEpipolarLineCoeff = F*leftPoint;
    a = rightEpipolarLineCoeff(1); 
    b = rightEpipolarLineCoeff(2); 
    c = rightEpipolarLineCoeff(3);
    
    % left epiline 
    xRight = [0; size(right,2)]; yRight = (-c-a*xRight) / b;
    
    % plot epipolar line in the right image
    plot(xRight, yRight, 'b-');
end

% show epipoles in a x-y plane
figure('name', 'Epipole e in the left image in red and Epipole e prime in the right image in blue');
plot(leftEpipolePosition(1),leftEpipolePosition(2),'ro');
hold on
plot(rightEpipolePosition(1),rightEpipolePosition(2),'bo');
xlabel('x');
ylabel('y');

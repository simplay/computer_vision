clear all
close all
clc

shouldRunSamples = true;
useLeft = false;

if shouldRunSamples
    if useLeft
        disp('running on given rotation Rl and translation Tl');
    else
        disp('running on given rotation Rr and translation Tr');
    end
else
    disp('running on estimated rotation R and translation t');
end

% Extrinsic parameters
Rr = [ 0.92848,  -0.12930,   0.34815;
       0.00000,   0.93744,   0.34815;
      -0.37139,  -0.32325,   0.87039];

Tr = [2.0;
      2.0;
      5.0];


Rl = [1.00000,   0.00000,   0.00000;
      0.00000,   0.92848,   0.37139;
      0.00000,  -0.37139,   0.92848];

Tl = [0.0;
      2.0;
      5.0];

left = mean(double(imread('Matched Points/left.jpg')),3);
right = mean(double(imread('Matched Points/Right.jpg')),3);

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

A=load('Matched Points/Matched_Points.txt');
[M, N] = size(A);

leftPoints = [A(:,3)'; A(:,4)'; ones(1,M)];
rightPoints = [A(:,1)'; A(:,2)'; ones(1,M)];

% estimate the fundamental matrix from the given matching points.
F = eightPointsAlgorithm(leftPoints,rightPoints); 
disp('Fundamental matrix F is equal to:')
disp(F);
if (rank(F) == 2)
    disp('F has rank 2');
else
    disp('Error: F has not rank 2');
end

% TODO: Compute Essential Matrix (Question 1)
K = [-83.33333,     0.00000,   250.00000;
       0.00000,   -83.33333,   250.00000;
       0.00000,     0.00000,     1.00000];

% computation of the essential matrix.
E = K'*F*K;
disp('Essential matrix E is equal to:')
disp(E);

% TODO: Compute Rotations and translatiosn between views (Question 2)

R90 = [0 -1 0; 1 0 0; 0 0 1];

[U,S,V] = svd(E);

%%

% compute all 4 candidate rotations
tmpR90 = R90;
sign = 1;
idx = 0;
candidateRotation = [];
for k=1:4
    tmpR = sign*U*tmpR90*V';
    
    % if determinate is greater zero.
    % actually for any rotation matrix, the determinate is either 1 or -1.
    if det(tmpR) > 0
        idx = idx + 1;
        candidateRotation = [candidateRotation,tmpR];
    end
    
    % sign swapping after first half of iterations
    if (k==2) sign = -sign; end
    
    % transpose matrix: every even iteration it is flipped initial version. 
    tmpR90 = tmpR90';
end

% reshape to a 3x3xN tensor condaining all N 3x3 candidate rotational transformations. 
candidateRotation = reshape(shiftdim(candidateRotation,1)', 3, 3, idx);

% set candidate translations
candidateTranslations = zeros(3,1,2);
candidateTranslations(:,:,1) = U(:,end);
candidateTranslations(:,:,2) = -U(:,end);

% compute height values from rotation and translation transformation
% combinations.
f = 4;
zComponents = zeros(M,1, size(candidateTranslations,3)*size(candidateRotation,3));
idx = 1;
for k=1:size(candidateRotation,3)
    for j=1:size(candidateTranslations,3)
        
        % override estimated transformations by given sample data.
        if shouldRunSamples
            if useLeft
                candidateRotation(:,:,k) = Rl;
                candidateTranslations(:,j) = Tl;
            else
                candidateRotation(:,:,k) = Rr;
                candidateTranslations(:,j) = Tr;
            end
        end
        
        % Apply the formula from en.wikipedia.org/wiki/Essential_matrix
        % in order to reconstruct the height z.
        r1 = repmat(candidateRotation(1,:,k),size(rightPoints,2),1);
        pr_x = repmat(rightPoints(1,:)',1,3);
        r3 = repmat(candidateRotation(3,:,k),size(rightPoints,2),1);
        r1_minus_prx_r3 = r1 - pr_x.*r3;
        t = repmat(candidateTranslations(:,j), 1, size(rightPoints,2))';
        z = sum(r1_minus_prx_r3.*t, 2)./sum(r1_minus_prx_r3 .* leftPoints',2);
        
        % save current height values.
        zComponents(:,:,idx) = z;
        idx = idx + 1;
    end
end

% take height values with most points in front of the camera,
% this directly relates to the best rotation and best translation.
bestZIdx = -1;
bestNumPosZ = -1;
for k=1:size(zComponents,3)
    currentZ = zComponents(:,:,k);
    
    numPositiveZ = length(find(currentZ >= 0));
    
    % relax values
    if(numPositiveZ > bestNumPosZ)
        bestNumPosZ = numPositiveZ;
        bestZIdx = k;
    end
end

%plot 3D points
z = zComponents(:,:,bestZIdx);
x = (z .* leftPoints(1,:)')/f;
y = (z .* leftPoints(2,:)')/f;

figure('name', 'Reconstructed 3D-image');
scatter3(x,y,z, 'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', [0.5,0.5,0.5]);

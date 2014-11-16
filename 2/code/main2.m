clear all
%close all
clc

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
[M N] = size(A);



leftPoints = [A(:,3)'; A(:,4)'; ones(1,M)];
rightPoints = [A(:,1)'; A(:,2)'; ones(1,M)];

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

E = K'*F*K;
disp('Essential matrix E is equal to:')
disp(E);




Rl = zeros(3);
tl = zeros(3,1);
Rr = zeros(3);
tr = zeros(3,1);
% TODO: Compute Rotations and translatiosn between views (Question 2)

R90 = [0 -1 0; 1 0 0; 0 0 1];


[U,S,V] = svd(E);
tmpR90 = R90;
sign = 1;
idx = 0;
candidateRotation = [];
%%

% compute all 4 candidate rotations
for k=1:4
    tmpR = sign*U*tmpR90*V;
    
    if det(tmpR) > 0
        idx = idx + 1;
        candidateRotation = [candidateRotation,tmpR];
    end
    
    if (k==2) sign = -sign; end
    tmpR90 = tmpR90';
end
candidateRotation = reshape(shiftdim(candidateRotation,1)', 3, 3, idx);

% TODO: Reconstrct the 3D points (Question 3)

% TODO: Visualize 3D points (Question 4)

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

% reshape to a 3x3xN tensor condaining all N 3x3 candidate rotational transformations. 
candidateRotation = reshape(shiftdim(candidateRotation,1)', 3, 3, idx);
R = candidateRotation;

candidateTranslations = zeros(3,1,2);
t_x = V*R90*S*V';
t = [t_x(3,2); t_x(1,3); t_x(2,1)];
t = t/norm(t);
candidateTranslations(:,:,1) = U(:,end);
candidateTranslations(:,:,2) = -U(:,end);
t = candidateTranslations;
% get best combination
current = -Inf;
f = 4;
for k=1:2
    for j=1:2
        %R(:,:,k) = Rr;
        %t(:,j) = Tr;
        a = repmat(R(1,:,k),size(rightPoints,2),1) - repmat(rightPoints(1,:)',1,3).*repmat(R(3,:,k),size(rightPoints,2),1);
        x3 = sum(a .* repmat(t(:,j), 1, size(rightPoints,2))',2) ./  sum(a .* leftPoints',2);
        
        % update points if we find a rotation, translation combination
        % where numbers in front of the camera is larger
        numPositiveZ = length(find(x3 >= 0));
        if numPositiveZ > current
            current = numPositiveZ;
            x1 = x3 ./ f .* leftPoints(1,:)';
            x2 = x3 ./ f .* leftPoints(2,:)';
            p = [x1,x2,x3]';
        end
    end
end


%plot 3D points
x = p(1,:)';
y = p(2,:)';
z = p(3,:)';

figure(2);
scatter3(x,y,z, 'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', [0.5,0.5,0.5]);
title('Reconstructed 3D-image');


%candidateTranslations = size(

% TODO: Reconstrct the 3D points (Question 3)

% TODO: Visualize 3D points (Question 4)

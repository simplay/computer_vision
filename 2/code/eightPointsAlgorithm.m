function F = eightPointsAlgorithm(x1, x2)
% function Fundamental_Matrix =Eight_Point_Algorithm(x1,x2)
% Calculates the Fundamental matrix between two views from the normalized 8 point algorithm
% Inputs: 
%               x1      3xN     homogeneous coordinates of matched points in view 1
%               x2      3xN     homogeneous coordinates of matched points in view 2
% Outputs:
%               F       3x3     Fundamental matrix

    [x1_bar, T1] = normalizePosition(x1);
    [x2_bar, T2] = normalizePosition(x2);
    
    
end


function [normalizedPositions, Transformation] = normalizePosition(positions)
    pixelArea = 2;
    N = length(positions);
    
    % since all positions are in homogenous coordinates
    % select only the first two rows
    xyPositions = positions(1:2,:);
    
    % compute average x and y coordinate which is supposed to define the 
    % center position of the given point cloud.
    avgPosition = (1/N)*sum(xyPositions, 2);
    
    % shift all positions towards the average position
    centeredXY = xyPositions - repmat(avgPosition, 1, N);
    
    % compute distance from all position to the average 'center' position.
    dist2Center = sqrt(centeredXY(1, :).^2 + centeredXY(2, :).^2);
    
    % compute the average distance to the average position among the
    % points in the point cloud.
    avgDistance = (1/N)*sum(dist2Center);
    
    % average distance to the points = sqrt(2)
    s = sqrt(pixelArea)/avgDistance;
    
    % scale all coordiantes by s
    % and shift them by scaled center
    Transformation = [s, 0, -s*avgPosition(1); 
         0, s, -s*avgPosition(2);
         0, 0, 1];
    
    % [(s*xy - repmat(s*avgPosition, 1, N));ones(1,N)]
    normalizedPositions = Transformation*positions;
end
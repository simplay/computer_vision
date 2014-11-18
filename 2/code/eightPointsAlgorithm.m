% CV-Project2
% Michael Single
% 08-917-445
function F = eightPointsAlgorithm(leftPositions, rightPositions)
% function Fundamental_Matrix =Eight_Point_Algorithm(x1,x2)
% Calculates the Fundamental matrix between two views from the normalized 8 point algorithm
% Inputs: 
%               x1      3xN     homogeneous coordinates of matched points in view 1
%               x2      3xN     homogeneous coordinates of matched points in view 2
% Outputs:
%               F       3x3     Fundamental matrix
    
    % Compute transformation matrices for left and right image feature positions.
    leftTransform = getPositionNormalizationTransform(leftPositions);
    rightTransform = getPositionNormalizationTransform(rightPositions);
    
    % transform points from left and right image to their mean position
    % having a distance of sqrt(2) to that position.
    normalizedLeftPositions = leftTransform*leftPositions;
    normalizedRightPositions = rightTransform*rightPositions;
    
    % define helper vectors in order to compute Y
    u = normalizedLeftPositions(1,:);
    v = normalizedLeftPositions(2,:);
    u_ = normalizedRightPositions(1,:);
    v_ = normalizedRightPositions(2,:);
    
    % Y is a homogenous matrix formed by all product combinations
    % of the tranformed (x,y) coordinates of both image positions.
    Y = [u_.*u; u_.*v; u_; v_.*u; v_.*v; v_; u; v];
    Y = [Y; ones(1,length(leftPositions))];
    
    % compute the fundamental matrix that is formed by Y's U rotation
    % matrix (resulting from Y's svd decomposition).
    [U,~,~] = svd(Y);
    F = reshape(U(:,end),3,3)';
    
    % ensure rank-2 constraint by setting the 3rd singular value to zero.
    [U, S, V] = svd(F);
    S(3,3) = 0;
    
    % reconstruct F (which is now rank two). Also normalze it by the image
    % transformations.
    F = U*S*V';
    
    % transform the fundamental matrix back to original units.
    F = rightTransform'*F*leftTransform;
end


function normalizationTransformation = getPositionNormalizationTransform(positions)
% Compute a Transformation Matrix such that when the given position inputs
% are transformed by this matrix they are approximately centered around the
% origin. Furthermore such points will have an average distance of sqrt(2)
% to the origin. Note that the origin is determined as the mean position of
% a given point cloud.
% @param position a 3 x N Matrix, containing N homogenous 2d positions.
% @return normalizationTransformation homohenous transform that transforms
%         points from the input to the origin (mean position in
%         pointcloud having an average distance of sqrt(2) to the origin).

    
    % get number of points
    N = length(positions);
    
    % since all positions are in homogenous coordinates
    % select only the first two rows - i.e. the (x,y) components of the points.
    xyPositions = positions(1:2,:);
    
    % compute mean x and y coordinate which is supposed to define the 
    % center position of the given point cloud.
    meanPosition = mean(xyPositions, 2);
    
    % shift all positions towards the average position, i.e. shift them to
    % the origin.
    centeredXY = xyPositions - repmat(meanPosition, 1, N);
    
    % compute the distance from all position to the mean 'center' position.
    dist2Center = sqrt(centeredXY(1, :).^2 + centeredXY(2, :).^2);
    
    % compute the mean distance to the mean position among the
    % points in the point cloud.
    meanDistance = mean(dist2Center);
    
    % ensure that the mean distance of the points to the origin is sqrt(2).
    % for this we define the scale factor s. I.e. center coordinates such
    % that mean squared distance between the origin and the points is 2
    % pixels.
    s = sqrt(2)/meanDistance;
    
    % define a homogenous transformation matrix that scales that
    % tranforms the input points to the origin (to their mean position)
    % having a distance of sqrt(2) to this origin.
    % scale all coordiantes by s and shift them by scaled center
    normalizationTransformation = [s, 0, -s*meanPosition(1); 
         0, s, -s*meanPosition(2);
         0, 0, 1];
end
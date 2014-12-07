function matches = getMatchingFeaturePoints(leftDescriptors, rightDescriptors)
%GETMATCHINGFEATUREPOINTS Summary of this function goes here
%   @param leftDescriptors  (N1 x 128) feature descriptor. 
%   @param rightDescriptors (N2 x 128) feature descriptor.
    
    % distances from each point to any other point
    squaredPointSetDistances = dist2(leftDescriptors, rightDescriptors);
    
    [smallesVals, indices] = min(squaredPointSetDistances, [], 2);
    for k=1:size(leftDescriptors, 1), 
        squaredPointSetDistances(k, indices(k)) = Inf; 
    end
    secSmallesVals = min(squaredPointSetDistances, [], 2);
    
    matchingIndices = find(1.5*smallesVals - secSmallesVals < 0);
    matches = [matchingIndices,indices(matchingIndices)];

end


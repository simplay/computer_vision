function matches = getMatchingFeaturePoints(leftDescriptors, rightDescriptors)
%GETMATCHINGFEATUREPOINTS Summary of this function goes here
%   @param leftDescriptors  (N1 x 128) feature descriptor. 
%   @param rightDescriptors (N2 x 128) feature descriptor.
    
    % distances from each point to any other point
    squaredPointSetDistances = dist2(leftDescriptors, rightDescriptors);
    
    % extract the current minima from the squared distance pairs
    % then again search for the minimum distance.
    [smallesVals, indices] = min(squaredPointSetDistances, [], 2);
    for k=1:size(leftDescriptors, 1), 
        squaredPointSetDistances(k, indices(k)) = Inf; 
    end
    secSmallesVals = min(squaredPointSetDistances, [], 2);
    
    % check whether the delta between the scaled minima and the other
    % values (represented by the 2nd smalles values). is smaller than one,
    % then we have found a strong minimum value, otherwise the determined
    % minimum is not signigicantly smaller than the other corresponding
    % values.
    % scaling:makes matches more robust - only take significant minima
    % according to a given threshold value which is here 1.
    matchingIndices = find(2.0*smallesVals - secSmallesVals < 0);
    matches = [matchingIndices,indices(matchingIndices)];

end


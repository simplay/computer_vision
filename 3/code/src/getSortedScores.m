function [scores, idxs] = getSortedScores( histSelFeatures, histData, frameCount)
%GETSORTEDSCORES Compute the similarity score between a given histogram and
%a set of other histograms. The similarity is determined by the normalized
%dot product.
%   @param histSelFeatures
%   @param histData
%   @param frameCount
    
    % compute dote product between selected histogram and any other
    % histogram. Note that this is unnormalized.
    scores = dot(repmat(histSelFeatures, frameCount, 1), histData, 2);
    
    % compute normaliation factor for norming the range of the scores.
    normalizationFScale = sqrt(sum(histData.^2, 2));
    normalFrameF = (repmat(norm(histSelFeatures), frameCount, 1).*normalizationFScale);
    
    % normalize the score and sort it according to its values in a
    % descending order.
    [scores, idxs] = sort(scores./normalFrameF, 'descend');
end


function [ histData ] = computeFrameHisto(descrCount, membership, frameCount, wordsPerCluster)
%COMPUTEFRAMEHISTO compute the histogram of how many feature poinst match a
%given cluster (identifier). i.e. it counts the matching features per
%cluster (i.e. word).
%   Detailed explanation goes her
    histData = zeros(frameCount, wordsPerCluster);
    memIdx = 0;
    for k=1:frameCount
        memEnd = memIdx + descrCount(k);
        memFrame = membership(memIdx+1:memEnd);
        memIdx = memEnd;
        histData(k,:) = histc(memFrame, 1:wordsPerCluster);    
    end

end


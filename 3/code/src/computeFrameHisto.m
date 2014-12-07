function [ histData ] = computeFrameHisto(descrCount, membership, frameCount, wordsPerCluster)
%COMPUTEFRAMEHISTO Summary of this function goes here
%   Detailed explanation goes here
    histData = zeros(frameCount, wordsPerCluster);
    memIdx = 0;
    for k=1:frameCount
        memEnd = memIdx + descrCount(k);
        memFrame = membership(memIdx+1:memEnd);
        memIdx = memEnd;
        histData(k,:) = histc(memFrame, 1:wordsPerCluster);    
    end

end


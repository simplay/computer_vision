function [ currentFrameName, featureCount, positions, orients, scales, descriptors ] = loadOwnSiftDataOf( siftMatFileName )
%LOADOWNSIFTDATAOF Summary of this function goes here
%   Detailed explanation goes here
    siftFilePathName = strcat('data/', siftMatFileName);
    
    load(siftFilePathName, 'currentFrameName', 'featureCount', ...
         'positions', 'orients', 'scales', 'descriptors'); 
end


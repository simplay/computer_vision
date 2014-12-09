function [ currentFrameName, featureCount, positions, orients, scales, descriptors ] = loadOwnSiftDataOf( siftMatFileName )
%LOADOWNSIFTDATAOF load sift data of given sift file
%   helper to parse the file name easily.
    siftFilePathName = strcat('data/', siftMatFileName);
    
    load(siftFilePathName, 'currentFrameName', 'featureCount', ...
         'positions', 'orients', 'scales', 'descriptors'); 
end


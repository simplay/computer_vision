function computeSiftDataOf(baseFrameName)
%COMPUTESIFTDATAOF generates the sift mat files containg all relevant sift
%data infromation for every frame of the passed base name. assembles also
%all the sift mat files int a one big _all.mat file.
%   @param String base name of frames extracted from provided video.
%          corresponds to the video name used.

    
    frameFiles = strcat('frames/',baseFrameName,'_*.png');
    frameFileNames = dir(frameFiles);
    
    fileCount = length(frameFileNames);
    allDescr = [];
    allPos = [];
    allScales = [];
    allOrients = [];
    allFrameNames = [];
    descrCount = zeros(fileCount, 1);
    
    
    for k = 1:fileCount
        currentFrameName = frameFileNames(k).name;
        currentGrayFrame = rgb2gray(imread(strcat('frames/',currentFrameName)));
        
        % currentGrayFrame is supposed to be a gray-scale image in single precision
        % Each column of 'featureFrames' is a feature frame
        % Each column of 'descriptors' is the descriptor of the corresponding frame
        [featureFrames, descriptors] = vl_sift(single(currentGrayFrame));
        
        featureCount = size(featureFrames, 2);
        positions = featureFrames(1:2,:)';
        orients = featureFrames(4,:)';
        scales = featureFrames(3,:)';
        descriptors = double(descriptors')/255;
        descrCount(k) = size(descriptors, 1);
        allDescr = [allDescr; descriptors];
        allPos = [allPos; positions];
        allOrients = [allOrients; orients];
        allScales = [allScales; scales];
        nameReps = size(descriptors,1);
        allFrameNames = [allFrameNames; repmat({currentFrameName}, nameReps, 1)];
        
        tillIdx = regexpi(currentFrameName, '\.');
        frameBaseName = currentFrameName(1:tillIdx-1);
        
        saveFilePathName = strcat('data/',frameBaseName, '.mat');
        save(saveFilePathName, 'currentFrameName', 'featureCount', ...
             'positions', 'orients', 'scales', 'descriptors'); 
        clc
        disp(['SIFT data  set for frame ', num2str(k), '/',num2str(fileCount), ' generated']); 
    end
    
    saveFilePathName = strcat('data/',baseFrameName, '_all.mat');
    save(saveFilePathName, 'allFrameNames', 'descrCount', 'allPos', ...
         'allOrients', 'allScales', 'allDescr');
    disp([num2str(fileCount), ' SIFT data sets from frames generated.']);
end


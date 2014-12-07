function extractVideoFrames(videoName)
%EXTRACTVIDEOFRAMES extracts all frames from a given video, stored in the
%folder 'video' and saves the frames as images in the folder 'frame'. as
%'videoName#frame.png'.
%   @param videoName String filename of video

    filePathName = strcat('video/', videoName);
    video = VideoReader(filePathName);
    
    % extract file name of video - without file extension.
    tillIdx = regexpi(videoName, '\.');
    baseFrameName = videoName(1:tillIdx-1);
    
    % extract frames as images from given video file.
    framePath = 'frames/';    
    for frameIdx=1:video.NumberOfFrames,
        currentFrame = read(video, frameIdx);
        currentFrameName = strcat(baseFrameName, '_', num2str(frameIdx), '.png');
        imwrite(currentFrame, strcat(framePath,currentFrameName));
    end

end


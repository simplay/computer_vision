function [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % Input:
  %  chromeDir (string) -- directory containing chrome images.
  %  nDir -- number of different light source images.
  %  chatty -- true to show results images. 
  % Return:
  %  L is a 3 x nDir image of light source directions.

  % Since we are looking down the z-axis, the direction
  % of the light source from the surface should have
  % a negative z-component, i.e., the light sources
  % are behind the camera.
    
  if ~exist('chatty', 'var')
    chatty = false;
  end
    
  mask = imread([chromeDir, 'chrome.mask.png']);
  mask = mask(:,:,1) / 255.0;

  for n=1:nDir
    fname = [chromeDir,'chrome.',num2str(n-1),'.png'];
    im = imread(fname);
    imData(:,:,n) = im(:,:,1);           % red channel
  end
  
  % number of rows is height, num. of columns is width.
  [imgHeight, imgWidth, imgCount] = size(imData);
  
  % sphere radius
  r = maskSphereRadius(mask, imgHeight, imgWidth);
  

  
  L = 1;
end

function r = maskSphereRadius(mask, imgHeight, imgWidth)
    % find (min,max)x(height,width) indices which are equal to one on mask
    % i.e. the boundary indices. compute their diameters and then, from
    % those the corresponding radius.
    
    % make sure min will be updated 
    minHeightIdx = imgHeight+ 1;
    minWidthIdx = imgWidth + 1;
    maxHeightIdx = -1;
    maxWidthIdx = -1;
    
    % find '1' boundaries
    for m=1:imgHeight,
        for n=1:imgWidth,
            if mask(m,n) == 1
                minHeightIdx = relax(m, minHeightIdx, @isSmaller);
                minWidthIdx = relax(n, minWidthIdx, @isSmaller);
                maxHeightIdx = relax(m, maxHeightIdx, @isGrater);
                maxWidthIdx = relax(n, maxWidthIdx, @isGrater);
            end
        end
    end
    % we have two diamenters: from minWidth to maxWidth
    % and from minHeight to maxHeight. 
    % We compute for these two diameters their radius - r = d/2 - 
    % and compute the mean radius between these two radi.
    diameters = [maxHeightIdx, maxWidthIdx]-[minHeightIdx, minWidthIdx];
    r = mean(diameters/2);
    
    % internal anonymous functions for having a simpler handler handling.
    function state = isSmaller(x,y)
        % tests x < y
        state = x < y;
    end

    function state = isGrater(x,y)
        % tests x > y
        state = x > y;
    end
end

function newOptimum = relax(candidateOptimum, currentOptimum, relation)
% if candidate optimum is better than current optimum according to a given
% relation, then update the current optimum by the candidate optimum.
% @param currentOptimum
% @param candidateOptimimum
% @param relation function handle tha returns a 1 for true & 0 for false.
    newOptimum = currentOptimum;
    if(relation(candidateOptimum, currentOptimum) == 1)
        newOptimum = candidateOptimum;
    end
end

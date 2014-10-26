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
  
  % sphere mean radius
  r = sphereRadiusOf(mask);
  
  % initialize light matrix and cameray ray.
  L = zeros(3, imgCount);
  
  % since we use orthographic projection
  % the direction of the camera from any point is (0,0,-1.)
  cameraRay = [0, 0, -1];
  
  % compute lights and sphere center of mask
  sphereCenter = sphereCenterOf(mask, ones(imgHeight, imgWidth), 0);
  for idx=1:imgCount
      currentImg = imData(:,:,idx);
      
      % center of bright spot - specular reflection - in current imgage.
      spotCenter = sphereCenterOf(currentImg, mask);
      
      % shift spot center to zero
      deltaCenter = spotCenter-sphereCenter;
      
      % compute sphere normal
      normal = computeSphereNormal(deltaCenter, r);
      
      % cosinus of incident angle: cosine-law: 
      % cos(theta) = dotProduct(sphere_normal, camera_ray)
      % note that sphere_normal, camera_ray have to be normalized.
      cosTheta = dot(normal, cameraRay);
      
      % According to Snell's reflection law:
      % given a position a camera ray hits the surface
      % in our case: the sphare surface, use the normal at this point.
      % refDir = camRay - 2*dot(camRay, normal)*normal
      L(:,idx) = 2*cosTheta*normal-cameraRay;
  end
end

function normal = computeSphereNormal(center, radius)
    % based on implicit equation of a sphere, i.e.
    % r^2 = x^2 + y^2 + z^2 
    % => z = sqrt(r^2 -(y^2 + z^2))
    % center contains shifted x,y coordinates.
    z = sqrt(radius^2 - sum(center.^2));
    
    % camera system convention
    normal = [center, -z];
    
    % normalization
    normal = normal/radius;
end


function center = sphereCenterOf(img, mask, threshold)
  % find center of specular reflection in a given image. No color-images.
  % assumption: there is only one specular reflection in an image.
  % @param img a mxn image containing a bright spot.
  % @param mask matrix for masking the given input image
  % @param threshold a positive integer, defining lower bound 
  %        for image pixel values defining a specular contribution.
  %        In other words: every pixel-value which has a value greater than
  %        this threshold can be considered as a specular contributing
  %        pixel. If this value is not passed, use a default threshold,
  %        relatively high - specular. Threshold equals zero are supposed
  %        to correspond to finding the center of a spherical mask.
  % @return coordinates in image (in pixel-coordinates) of specular center.
  if nargin < 3
    threshold = 250;  
  end
  
  [height, width] = size(img);
  
  % specular reflection - brightness spot
  % are assumed to be values bigger than our brigthness threshold
  acceptedValuesMask = (img-threshold) > 0.0;
  maskedImg = double(img).*double(double(mask).*acceptedValuesMask);
  
  % normalization factor: divide final result by all elements 
  % which were selected by using our mask(s). 
  nF = sum(maskedImg(:));
  
  % image indices tensors (colum,row)-indices.
  yIdx = repmat(1:height, width, 1)';
  xIdx = repmat(1:width, height, 1);
  
  % select all required indices determined by our masks.
  maskedXidx = maskedImg.*xIdx;
  maskedYidx = maskedImg.*yIdx;
  
  % sample center
  center = [sum(maskedXidx(:)), sum(maskedYidx(:))] / nF;
end

function r = sphereRadiusOf(mask)
% @param mask (mxn) boolean (i.e. 0-1-valued) matrix
    % find (min,max)x(height,width) indices which are equal to one on mask
    % i.e. the boundary indices. compute their diameters and then, from
    % those the corresponding radius.
    
    % get resolution of mask matrix
    [height, width] = size(mask);
    
    % make sure min will be updated 
    minHeightIdx = height+ 1;
    minWidthIdx = width + 1;
    maxHeightIdx = -1;
    maxWidthIdx = -1;
    
    % find '1' boundaries in mask
    for m=1:height,
        for n=1:width,
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

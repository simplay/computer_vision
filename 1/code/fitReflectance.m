function [n, albedo] = fitReflectance(img, L)
  % [n, albedo] = fitReflectance(im, L)
  % 
  % Input:
  %   img - nPix x nDirChrome array of brightnesses,
  %   L  - 3 x nDirChrome array of light source directions.
  % Output:
  %   n - nPix x 3 array of surface normals, with n(k,1:3) = (nx, ny, nz)
  %       at the k-th pixel.
  %   albedo - nPix x 1 array of estimated albdedos
    
  % for every pixel p in image img it holds true:
  % img(x) = albedo(x)*(normal_at(x)*L)
  
   %% [imgHeight, imgWidth , imgCount] = size(img);
   S = L';
   I = img';
   n_bar = (S'*S) \ (S'*I);
   n_bar = n_bar';
   
   albedo = sqrt(sum(n_bar.^2, 2));
  
   n = n_bar ./ repmat((albedo), 1, 3);

end



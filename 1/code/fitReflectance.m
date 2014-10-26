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
   
   % use same convention as introduced during the lecture.
   S = L'; I = img';
   
   % system using the pseudo inverse.
   n_tilde = (S'*S) \ (S'*I);
   n_tilde = n_tilde';
   
   % albedo is norm of 
   albedo = sqrt(sum(n_tilde.^2, 2));
   
   % normalize by dividing by albedo 
   n = n_tilde ./ repmat(albedo, 1, 3);
end
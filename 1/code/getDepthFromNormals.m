function [depth] = getDepthFromNormals(n, mask)
  % [depth] = getDepthFromNormals(n, mask)
  %
  % Input:
  %    n is an [N, M, 3] matrix of surface normals (or zeros
  %      for no available normal).
  %    mask logical [N,M] matrix which is true for pixels
  %      at which the object is present.
  % Output
  %    depth an [N,M] matrix providing depths which are
  %          orthogonal to the normals n (in the least
  %          squares sense).
  %
  [x_max, y_max, ~] = size(mask);
  v = []; i = []; j = []; s = [];
  
  % collect all indices for constraints system of eqn.
  row_idx = 1;
  for y=1:y_max
      for x=1:x_max
          if mask(x,y) == 1
              
              % tangent t_x = (1,0,Z_x)'
              % Z_x ~= Z(x+1, y) - Z(x,y)
              j = [j,row_idx, row_idx];
              i = [i, rowIndexOf(x + 1, y, x_max), rowIndexOf(x, y, x_max)];
              s = [s, n(x,y,3), - n(x,y,3)];
              row_idx = row_idx + 1;
              
              % tangent t_y = (0,1,Z_y)'
              % Z_y ~= Z(x, y+1) - Z(x,y)         
              j = [j,row_idx, row_idx];
              i = [i, rowIndexOf(x, y + 1, x_max), rowIndexOf(x, y, x_max)];
              s = [s, n(x,y,3), - n(x,y,3)];
              row_idx = row_idx + 1;
              
              % inverts x,y direction
              v = [v; -n(x,y,2); -n(x,y,1)];
           end
      end
  end
  
  % set pixel at top left to zero depth
  j = [j, row_idx];
  i = [i, rowIndexOf(1, 1, x_max)];
  s = [s, 1];
  v = [v; 0];
  A = sparse(j,i,s, max(j), x_max*y_max);
  
  % solve for z of System Az = v
  z = A\v;
  depth = reshape(z, [x_max, y_max, 1]);
end

function index = rowIndexOf(x, y, max)
    index = max*(y - 1) + x;
end

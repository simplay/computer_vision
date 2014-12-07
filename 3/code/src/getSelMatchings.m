function [ onPos, onOrients, onScales ] = getSelMatchings(matchIdxs, positions, orients, scales, maskedIndices)
%GETSELMATCHINGS Summary of this function goes here
%   Detailed explanation goes here
    
    % temporary assignment
    tPos = positions;
    tOrients = orients;
    tScales = scales;
    
    % are we dealing with the left (image selected region comes from).
    if(nargin > 4)
      tPos = positions(maskedIndices,:);
      tOrients = orients(maskedIndices,:);
      tScales = scales(maskedIndices,:);
    end
    
    % select matched feature attributes.
    onPos = tPos(matchIdxs, :);
    onOrients = tOrients(matchIdxs, :);
    onScales = tScales(matchIdxs, :);
end


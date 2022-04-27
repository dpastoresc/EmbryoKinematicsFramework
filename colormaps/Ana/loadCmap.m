% Load colormap
% cmap = loadCmap(name,levels)
%
% Available colormaps (specify 128 or 256):
% 
% '5_ramps'
% 'edges'
% 'ice'
% 'phase'
% 'redblue'
% 'split_bluered_warmmetal'
%

function [cmap, cmapor]= loadCmap(name,levels)

cmap_name = [name '-' num2str(levels) '.jpg'];
% Read colormap
[A, cmapor] = imread(cmap_name);
% Remove singleton dimension
A = squeeze(A);
% We have to rescale the values between 0 and 1
% Our color levels are from 0-255, independently of im_length
% scaledI = (I-min(I(:))) ./ (max(I(:)-min(I(:))));
cmap = double(A)./255;

% demo_colormap_signed.m
% 
% This script demonstrates use of colormap_signed. For demonstration
% purposes, it uses a modification of Matlab's 'peaks' surface consisting 
% of shifted, scaled Gaussians. Peaks was modified to produce a surface where
% the magnitudes of the maximum and minimum curvatures are not equal. In this
% case, applying a symmetrical signed colormap would assign the 'zero color'
% to a nonzero value of curvature (the value halfway between the min and the
% max). The colormap_signed function lets the user specify the 'zero
% location' that corresponds to the 'zero color' (i.e., white) on the colormap.
% The colormap is stretched/compressed to use the full color range for the data
% while ensuring the that the 'zero color' correctly labels regions of zero 
% value (curvature, for this example).
 
% by Peter Hammer, April 2015

[x,y] = meshgrid(linspace(-3,3,256));
p =  4*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ...
   - 3*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2);

k=del2(p,6/255); % discrete approximation of Laplacian

% compute relative location of zero value (of curvature) between min and
% max values of curvature
zero_loc=(0-min(k(:)))/(max(k(:))-min(k(:)));

% generate colormap (128 rows x 3 columns)
cmap=colormap_signed(128,zero_loc);

% plot surface and apply colormap
figure('color',[.7 .7 .7])
surf(p,k,'facecolor','interp',...
    'edgecolor','none',...
    'ambientstrength',0.4,...
    'specularstrength',0.2)
light('position',[0 0 10])
light('position',[0 0 -10])
axis off
colormap(cmap)
colorbar
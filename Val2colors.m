function Colors = Val2colors(txt,cl);
%
% Syntax :
% Colors = Val2colors(txt,cl);
%
% This function creates a colormap table according to the values of the txt
% variable. The colormap can be also specified.
%
% Input Parameters:
%      txt        : Vector of characteristic values
%       cl        : Colormap
%
%
% Output Parameter:
%     Colors      : Color Table
%
% See also: Surf_Color
%__________________________________________________
% Authors: Yasser Aleman Gomez
%LIM
% February 1th 2013
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ==1
    cl = 'jet';
end
[ut,i1t,it2] = unique(txt);
switch cl
    case 'hot'
        col = hot(size(ut,1));
    case 'hsv'
        col = hsv(size(ut,1));
    case 'jet'
        col = jet(size(ut,1));
    case 'cool'
        col = cool(size(ut,1));
    case 'bone'
        col = bone(size(ut,1));
    case 'pink'
        col = pink(size(ut,1));
    case 'winter'
        col = winter(size(ut,1));
    case 'autumn'
        col = autumn(size(ut,1));
    case 'summer'
        col = summer(size(ut,1));
    case 'spring'
        col = spring(size(ut,1));
    case 'copper'
        col = copper(size(ut,1));
    case 'spectral'
        col = spectral(size(ut,1));
    case 'gray'
        col = gray(size(ut,1));
end
Colors = col(it2,:);
return
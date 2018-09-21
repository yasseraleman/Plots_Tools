function Colors = Image2Colors(I,cl);
%
% Syntax :
% Colors = Image2Colors(I,cl);
%
% The following program converts the Image I into a colored image named
% Colors. Eack voxel in I is coloured acording its intensity value
% following the RGB protocol.
%
% Input Parameters:
%         I            : 3D Image (Matlab variable NxMxP)
%        cl            : Color map (ie. jet, hsv, hot, pink, etc);
%
%
% Output Parameters:
%        Colors        : Coloured Image
%
% Related references:
%
% See also:
%__________________________________________________
% Authors:  Yasser Aleman Gomez
% Unidad de Medicina y Cirugía Experimental, Hospital General Universitario Gregorio Marañón, Madrid
%
% November 1st 2012
% Version $1.0

ind = find(I==0);
txt = I(:);
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
r = reshape(Colors(:,1),[size(I)]);
r(ind) = 0;
g = reshape(Colors(:,2),[size(I)]);
g(ind) = 0;
b = reshape(Colors(:,3),[size(I)]);
b(ind) = 0;
if length(size(I))==2
    Colors = cat(3,r,g,b);
elseif length(size(I))==3
    Colors = cat(4,r,g,b);
end
return
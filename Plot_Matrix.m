function [hf] = Plot_Matrix(Mat,Nv,cl);
% Syntax :
% [Colors] = Plot_Matrix(Mat,Nv,cl);
%
% This function plots the connectivity matrix with white columns and rows dividing intraregional correlation matrices
%
% Input Parameters:
%   Mat         : Correlation Matrix
%   Nv          : Number of variables to create intra-region correlation matrix
%   cl          : Colormap used to see the results.
%
% Output Parameters:
%  hf           : Figure Handle.
%
% Related references:
% 
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% September 29th 2007
% Version $1.0

if nargin <3
    cl = 'jet';
    
end
if nargin <2
    errordlg('Wrong inputs');
end
%=========================Main program====================================%
A = zeros(size(Mat)+size(Mat)/Nv+1);
Nums = [1:Nv+1:size(A,1)];
cols = [1:size(A,1)];
cols(Nums) = [];
A(cols,cols) = Mat;

txt=A(:);
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
    case 'cyanyellow'
        col = cyan_yellow(size(ut,1));
end

[X,Y]= meshgrid(Nums,[Nums(1):Nums(end)]);X = X(:);Y = Y(:);
ind = sub2ind(size(A),[X;Y],[Y;X]);

Colors = col(it2,:);
R = zeros(size(A));
R(1:size(Colors,1)) = Colors(:,1);
    
G = zeros(size(A));
G(1:size(Colors,1)) = Colors(:,2);

B = zeros(size(A));
B(1:size(Colors,1)) = Colors(:,3);

R(ind) = 1;
B(ind) = 1;
G(ind) = 1;
C = cat(3,R,G,B);

hf = figure('numbertitle','off','Color','white','Position',[0 0 1200 900]);
imshow(C,'InitialMagnification','fit');
h = colorbar;
siz = get(h,'Ylim');
Yticks = siz(1):(siz(2) - siz(1))/5:siz(2);
set(h,'YTick', Yticks);
YLabs = linspace(min(txt),max(txt),6);
set(h,'YTickLabel', cellstr(num2str(YLabs(:))));

return;
function varargout = colormaps_colors(varargin);
% Syntax :
% [col] = colormaps_colors(cl,n);
%
% This function generates a colormaps. Tipying colormaps_colors without
% input in the command windows to displays available colormaps.
%
% Input Parameters:
%   cl          : Colormap name.
%   n           : Number of colors
%   cent_loc    : Location of central color. Important for signed
%                 colormaps.
%
% Output Parameters:
%   col         : Output colormap n x 3 matrix.
%
% Related references:
% 
%
% See also: Smooth_Surf Surf_Comp Plot_Surf Plot_oversurf Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% September 29th 2007
% Version $1.0
%% ============================= Checking Inputs ======================= %%
if nargin == 0

    varargout{1} = strvcat('jet',...
            'hot',...
            'bluered',...
            'greenred',...
            'hsv',...
            'cool',...
            'bone',...
            'pink',...
            'winter',...
            'autumn',...
            'summer',...
            'spring',...
            'copper',...
            'spectral',...
            'gray',...
            'cyanyellow',...
            'yellowcyan',...
            'green',...
            'red',...
            'blue',...
            'redblue',...
            'redgreen',...
            'belle',...
            'manitou',...
            'munising',...
            'saginaw',...
            'saugatuck',...
            'blackred',...
            'blackgreen',...
            'blackblue',...
            'blackyellow',...
            'blackcyan',...
            'blackmagenta',...
            'g-w-r',...
            'dkb-w-dkr',...
            'b-w-r',...
            'c-w-y');
        
%         disp(' ');
%         disp('============== Available Colormaps ================');
%         disp( varargout{1});
%         disp('===================================================');
        return;
elseif nargin == 1
    cl = varargin{1};
    n = 100;
    cent_loc = .5;
elseif nargin == 2
    cl = varargin{1};
    n = varargin{2};
    cent_loc = .5;
elseif nargin == 3
    cl = varargin{1};
    n = varargin{2};
    cent_loc = varargin{3};
else
    error('To Many Input Parameters');
    return;
end
if nargout > 1
    error('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%


%=========================Main program====================================%
switch cl
    case 'hot'
        col = hot(n);
    case 'hsv'
        col = hsv(n);
    case 'jet'
        col = jet(n);
    case 'cool'
        col = cool(n);
    case 'bone'
        col = bone(n);
    case 'pink'
        col = pink(n);
    case 'winter'
        col = winter(n);
    case 'autumn'
        col = autumn(n);
    case 'summer'
        col = summer(n);
    case 'spring'
        col = spring(n);
    case 'copper'
        col = copper(n);
    case 'spectral'
        col = spectral(n);
    case 'gray'
        col = gray(n);
    case 'cyanyellow'
        [col] = colorGradient([0 200 255],[255 200 0],n);
    case 'yellowcyan'
        [col] = colorGradient([255 200 0],[0 200 255],n);
    case 'green'
        [col] = colorGradient([0 30 0],[0 230 0],n);
    case 'red'
        [col] = colorGradient([30 0 0],[230 0 0],n);
    case 'blue'
        [col] = colorGradient([0 0 30],[0 0 230],n);
    case 'bluered'
        [col] = colorGradient([0 124 200],[210 50 26],n);
    case 'redblue'
        [col] = colorGradient([210 50 26],[0 124 200],n);
    case 'greenred'
        [col] =colorGradient([0 230 0],[230 0 0],n);
    case 'redgreen'
        [col] =colorGradient([230 0 0],[0 230 0],n);
    case 'belle'
        col = belle(n);
    case 'manitou'
        col = manitou(1,n);
    case 'munising'
        col = munising(1,n);
    case 'saginaw'
        col = saginaw(1,n);
    case 'saugatuck'
        col = saugatuck(1,n);
    case 'blackred'
        col = colorGradient([20 0 0],[255 0 0],n);
    case 'blackgreen'
        col = colorGradient([0 20 0],[0 255 0],n);
    case 'blackblue'
        col = colorGradient([0 0 20],[0 0 255],n);
    case 'blackyellow'
        col = colorGradient([20 20 0],[255 255 0],n);
    case 'blackcyan'
        col = colorGradient([0 20 20],[0 255 255],n);
    case 'custombone'
        col = colorGradient([156 180 180],[0 20 240],n);
   case 'blackmagenta'
        col = colorGradient([20 0 20],[255 0 255],n);
    case 'g-w-r'
        col = signed_colormap('g-w-r',n, cent_loc);
    case 'dkb-w-dkr'
        col = signed_colormap('dkb-w-dkr',n, cent_loc);
    case 'b-w-r'
        col = signed_colormap('b-w-r',n, cent_loc);
    case 'c-w-y'
        col = signed_colormap('c-w-y',n, cent_loc);
    otherwise
        col = jet(n);
end
%========================End of main program==============================%
% Outputs
varargout{1} = col;
return

function s = spectral(m)
%SPECTRAL Black-purple-blue-green-yellow-red-white color map.
%
%         map = spectral(num_colors)
%
% SPECTRAL(M) returns an M-by-3 matrix containing a "spectral" colormap.
% SPECTRAL, by itself, is the same length as the current colormap.
%
% For example, to reset the colormap of the current figure:
%
%           colormap(spectral)
%
% See also HSV, GRAY, PINK, HOT, COOL, BONE, COPPER, FLAG,
%          COLORMAP, RGBPLOT.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
base = [
    0.0000 0.0000 0.0000
    0.4667 0.0000 0.5333
    0.5333 0.0000 0.6000
    0.0000 0.0000 0.6667
    0.0000 0.0000 0.8667
    0.0000 0.4667 0.8667
    0.0000 0.6000 0.8667
    0.0000 0.6667 0.6667
    0.0000 0.6667 0.5333
    0.0000 0.6000 0.0000
    0.0000 0.7333 0.0000
    0.0000 0.8667 0.0000
    0.0000 1.0000 0.0000
    0.7333 1.0000 0.0000
    0.9333 0.9333 0.0000
    1.0000 0.8000 0.0000
    1.0000 0.6000 0.0000
    1.0000 0.0000 0.0000
    0.8667 0.0000 0.0000
    0.8000 0.0000 0.0000
    0.8000 0.8000 0.8000
    ];
n = length(base);
X0 = linspace (1, n, m);
s = interp1(1:n,base,X0);

return;
function col = signed_colormap(cl,n,cent_loc);



% % % % % % % % c = [0 1 1;...
% % % % % % % %     0 0 0.5;...
% % % % % % % %     0 0 1;...
% % % % % % % %     .9 .9 .9;...
% % % % % % % %     1 0 0;...
% % % % % % % %     0.5 0 0;...  
% % % % % % % %     1 1 0];





switch cl
    case 'c-w-y'
        c = [0 1 1;...
    0 0 0.5;...
    0 0 1;...
    .9 .9 .9;...
    1 0 0;...
    0.5 0 0;...  
    1 1 0];
    case 'b-w-r'
        c = [0 0 0.5;...
    0 0 1;...
    0 1 1;...
    1 1 1;...
    1 1 0;...
    1 0 0;...
    0.5 0 0];
    case 'dkb-w-dkr'
        c = [0 0 30;...
            0 51 102;...
            102 178 255;...
            240 240 240;...
            255 195 155;...
            237 97 0;...
            51 0 0]/255;
    case 'g-w-r' % Green-White-Red
        c = [0 1 0;...
            0 0.5 0;...
            0 0.25 0;...
            1 1 1;...
            .25 0 0;...
            0.5 0 0;...
            1 0 0];
end

i_mid = 0.5*(1+size(c,1));
cmap_neg=c(1:i_mid,:);
cmap_pos=c(i_mid:end,:);
i0 = 1+ round(n * cent_loc); % row of cmap (n rows) corresponding to zero 

x=(1:i_mid)'/i_mid;
cmap_neg_i=interp1(x,cmap_neg,linspace(x(1),1,i0));
cmap_pos_i=interp1(x,cmap_pos,linspace(x(1),1,n-i0));
col = [cmap_neg_i; cmap_pos_i];
return;
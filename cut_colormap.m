function varargout = cut_colormap(varargin);
% Syntax :
% s = cut_colormap(cl,n,N,boolflip);
%
% This function cuts a specified colormap.
%
% Input Parameters:
%   cl          : Colormap used to see the results.
%   n           : Number of points in the colormap
%   N           : Useful percent in the colormap
%   boolflip    : Boolean ariable to flip or not the colormap.
%
% Output Parameters:
%  Colors       : Output colormap matrix.
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

%=========================Main program====================================%
%% ============================= Checking Inputs ======================= %%
if nargin == 0
    cl = 'jet';
    n = 0;
elseif nargin == 1
    cl = varargin{1};
    n = 100;
    N = [0 100];
    boolflip = 0;
elseif nargin == 2
    cl = varargin{1};
    n = varargin{2};
    N = [0 100];
    boolflip = 0;
elseif nargin == 3
    cl = varargin{1};
    n = varargin{2};
    N = varargin{3};
    boolflip = 0;
elseif nargin == 4
    cl = varargin{1};
    n = varargin{2};
    N = varargin{3};
    boolflip = varargin{4};
elseif nargin > 4
    error('To Many Input Parameters');
    return;
end
if nargout > 1
    error('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%
%% ============================== Main Program ========================= %%
if length(N)< 2
    N = [0 N(:)];
end

[col] = colormaps_colors(cl,21);
Number = [floor(21*N(1)/100) floor(21*N(2)/100)];
if Number(1) == 0
    Number(1) = 1;
end
tempc = col(Number(1):Number(2),:);
if boolflip 
  tempc = flipdim(tempc,1);
end
base = [tempc];
m = length(base);
X0 = linspace (1, m, n);
s = interp1(1:m,base,X0);
varargout{1} = s;
%% ========================= End of Main Program ======================= %%
return
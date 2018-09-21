function cmap = colormap_signed(n,zero_loc)

% Construct colormap for displaying signed data. The function outputs an n x 3
% colormap designed for use with signed data. The user can specify the location
% in the data range that corresponds to zero, and the colormap is then constructed
% so that white maps to zero. 
%
% Input arguments:
%   n: number of rows in colormap (default = 64)
%   zero_loc: location of zero (fractional dist between neg and pos limits
%             of data). If k is a signed function of 2 variables that spans
%             a range from negative to positive, compute the location
%             of the zero value: 
%                   zero_loc = (0 - min(k(:)))/(max(k(:)) - min(k(:))) 
%
% Usage:
% cmap = colormap_signed returns a 64 x 3 colormap in which the middle
% rows tend toward white; lower rows (negative values) tend toward blue
% while higher rows tend toward red. Variables n and zero_loc assume default
% values of 64 and 0.5, respectively.
%
% cmap = colormap_signed(n) returns a n x 3 colormap, otherwise similar to
% above. Variable zero_loc assumes default value of 0.5. 
%
% cmap = colormap_signed(n,zero_loc) returns a n x 3 colormap in which the 
% location of the row corresponding to the 'zero color' is given by zero_loc.
% See section above on input arguments for example of how to compute zero_loc.
%
% NOTE: As the value of zero_loc deviates from 0.5, the colormap created by 
% this function becomes progressively warped so that the portion of the data 
% range mapped to warm colors does not equal that mapped to cool colors. This 
% is intentional and allows the full range of colors to be used for a given
% signed data range. If you want a signed colormap in which the incremental
% color change is constant across the entire data range, one option is to 
% use this function to return a symmetrical signed colormap (i.e., zero_loc 
% = 0.5) and then manually set the colorbar properties to crop the colorbar.

% Written by Peter Hammer, April 2015 and posted on Matlab File Exchange

switch nargin
    case 2
        if (n < 1)
            error('First input argument must be greater than zero.')
        end
        if ((zero_loc < 0) || (zero_loc > 1))
            error('Second input argument must be between 0 and 1.')
        end
    case 1
        zero_loc = 0.5;
        if (n < 1)
            error('First input argument must be greater than zero.')
        end
    case 0
        zero_loc = 0.5;
        n = 64;
    otherwise
        error('Too many input arguments.')
end

% Array c must have odd number of rows with 'zero color' in middle row.
% This is a modified jet colormap with white replacing green in the
% middle (DarkBlue-Blue-Cyan-White-Yellow-Red-DarkRed).
% % % % % % % % c = [0 0 0.5;...
% % % % % % % %     0 0 1;...
% % % % % % % %     0 1 1;...
% % % % % % % %     1 1 1;...
% % % % % % % %     1 1 0;...
% % % % % % % %     1 0 0;...
% % % % % % % %     0.5 0 0];

% % % % % % % % c = [0 1 1;...
% % % % % % % %     0 0 0.5;...
% % % % % % % %     0 0 1;...
% % % % % % % %     .9 .9 .9;...
% % % % % % % %     1 0 0;...
% % % % % % % %     0.5 0 0;...  
% % % % % % % %     1 1 0];





c = [0 0 51;...
    0 51 102;...
    102 178 255;...
    240 240 240;...
    255 195 155;...
    237 97 0;...  
    51 0 0]/255;

% Green-White-Red
% c = [0 1 0;...
%     0 0.5 0;...
%     0 0.25 0;...
%     1 1 1;...
%     .25 0 0;...
%     0.5 0 0;...  
%     1 0 0];







i_mid = 0.5*(1+size(c,1));
cmap_neg=c(1:i_mid,:);
cmap_pos=c(i_mid:end,:);
i0 = 1+ round(n * zero_loc); % row of cmap (n rows) corresponding to zero 

x=(1:i_mid)'/i_mid;
cmap_neg_i=interp1(x,cmap_neg,linspace(x(1),1,i0));
cmap_pos_i=interp1(x,cmap_pos,linspace(x(1),1,n-i0));
cmap = [cmap_neg_i; cmap_pos_i];
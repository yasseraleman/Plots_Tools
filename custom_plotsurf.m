function strsurf = custom_plotsurf(Surf, parameters);
%
% Syntax :
% strsurf = custom_plotsurf(Surf, parameters);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   Surfa       : Surface matrix.
%               : An atlas surface file is considered as a single surface.
%   parameters  : Patch parameters.
%
% Output Parameters:
%   strsurf     : Patch handle.
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Red_Surf Plot_oversurf Atlas_Surf
% Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 30th 2006
% Version $1.0

%=====================Checking Input Parameters===========================%
if nargin < 2
    parameters = {'EdgeAlpha', 1,...
        'EdgeColor','none',...
    'LineStyle', '-',...
    'LineWidth', .3,...
    'Marker', 'none',...
    'MarkerEdgeColor', 'auto',...
    'MarkerFaceColor', 'none',...
    'MarkerSize', 2,...
    'FaceAlpha',1,...
    'FaceLighting', 'gouraud',...
    'EdgeLighting', 'none',...
    'BackFaceLighting', 'reverselit',...
    'AmbientStrength', .3,...
    'DiffuseStrength', .6,...
    'SpecularStrength', .9,...
    'SpecularExponent', 60,...
    'SpecularColorReflectance',1,...
    'Visible','on'};
end
%=========================================================================%

%=========================Main Program====================================%
Surf.SurfData.FaceColor = 'interp';
if ~isfield(Surf.SurfData,'FaceVertexCData')
    strsurf = patch(Surf.SurfData,'facecolor',[.9 .9 .9],'edgecolor','none', 'tag','model0','facelighting','phong', parameters{:});
else
    strsurf = patch(Surf.SurfData,'edgecolor','none', 'tag','model0','facelighting','phong', parameters{:});
end
%=========================================================================%
return

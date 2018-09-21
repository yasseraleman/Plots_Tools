function Plot_Tracts_VTK(Tracts,config);
%
% Syntax :
% Plot_Tracts_VTK(Tracts,config);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   Surfa       : Tracks files (Tract according to trackvis format).
%   config      : Configuration Structure
%         % Plot_Surf Options   
%         config.tr   :      % Opacity vector
%         config.sa   :      % Plotting in the same window
%         config.cl   :      % Default Colormap   
%         
%         % MatVTK Options
%         config.ambient            :     % Ambient
%         config.ambientColor       :     % Ambient Color
%         config.backgroundColor    :     % Background color
%         config.cameraFocalDisk    :     % Camera Focal Disk
%         config.cameraFocalPoint   :     % Camera Focal Point
%         config.cameraPosition     :     % Camera Position
%         config.cameraViewUp       :     % Camera View Up
%         config.color              :     % Surface Color
%         config.diffuse            :     % Diffuse
%         config.diffuseColor       :     % Diffuse Color
%         config.edgeColor          :     % Edge Color
%         config.edgeVisibility     :     % Edge Visibility
%         config.interpolation      :     % Interpolation
%         config.lightFollowCamera  :     % Light Following the camera
%         config.lineWidth          :     % Line Width
%         config.opacity            :     % Opacity
%         config.orientation        :     % Orientation
%         config.origin             :     % Origin
%         config.position           :     % Position
%         config.representation     :     % Object Representation
%         config.scale              :     % Scale
%         config.shading            :     %  Shading
%         config.specular           : 0.6440;     % Specular percent of lighting
%         config.specularColor      : [0.80 0.80 0.80];    % Specular Color
%         config.specularPower      : 100;      % Specular power
%         config.twoSidedLighting   : 1;        % Two Side Lighting
% 
% Output Parameters:
%
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
if nargin==0
    [Tracts,sts] = spm_select([1 Inf],'any','Selecting Tracts Files','',cd);
end

if nargin <2
        % Tracts Options
        config.lineWidth          = 1;  % Sets the line width, defaults to 1.0
        config.lineTubes          = 1;  % 0=draw lines, 1=draw tubes. defaults to 0.
        config.tubeSides          = 8;  % Number of sides for the tube, defaults to 6 (hexagonal surface) 
        config.tubeRadius         = 0.5;  % Tube thickness, defaults to 1.
    
        % Plot_Surf Options   
        config.tr                 = 1;            % Opacity vector
        config.sa                 = 'y';          % Plotting in the same window
        config.cl                 = 'jet';        % Default Colormap   
        
        % MatVTK Options
        config.ambient            = 0;            % Ambient
        config.ambientColor       = [1 1 1];      % Ambient Color
        config.backgroundColor    = [0 0 0];      % Background color
        config.cameraFocalDisk    = 1;            % Camera Focal Disk
        config.cameraFocalPoint   = [0 0 0];      % Camera Focal Point
        config.cameraPosition     = [29.0283 26.0049 164.4381]; % Camera Position
        config.cameraViewUp       = [0 0 0];    % Camera View Up
        config.color              = [1 1 1];      % Surface Color
        config.diffuse            = 1;            % Diffuse
        config.diffuseColor       = [1 1 1];      % Diffuse Color
        config.edgeColor          = [0 0 0];      % Edge Color
        config.edgeVisibility     = 0;            % Edge Visibility
        config.interpolation      = 'gouroud';    % Interpolation
        config.lightFollowCamera  = 1;            % Light Following the camera
        config.lineWidth          = 1;            % Line Width
        config.opacity            = 1;            % Opacity
        config.orientation        = [0 0 0];      % Orientation
        config.origin             = [0 0 0];      % Origin
        config.position           = [0 0 0];      % Position
        config.representation     = 'surface';    % Object Representation
        config.scale              = [1 1 1];      % Scale
        config.shading            = 0;            %  Shading
        config.specular           = 0.6440;       % Specular percent of lighting
        config.specularColor      = [0.80 0.80 0.80];    % Specular Color
        config.specularPower      = 100;          % Specular power
        config.twoSidedLighting   = 1;            % Two Side Lighting
else

    % Tracts Options
    if ~isfield(config,'lineWidth')
        config.lineWidth          = 1;  % Sets the line width, defaults to 1.0
    end
    if ~isfield(config,'lineTubes')
        config.lineTubes          = 1;  % 0=draw lines, 1=draw tubes. defaults to 0.
    end
    if ~isfield(config,'tubeSides')
        config.tubeSides          = 4;  % Number of sides for the tube, defaults to 6 (hexagonal surface)
    end
    if ~isfield(config,'tubeRadius')
        config.tubeRadius         = 1;  % Tube thickness, defaults to 1.
    end
    
    % Plot_Surf configuration Parameters
    if ~isfield(config,'tr')
        config.tr                 = 1;             % Opacity vector
    end
    if ~isfield(config,'sa') 
        config.sa                 = 'y';           % Plotting in the same window
    end
    if ~isfield(config,'cl')
        config.cl                 = 'jet';         % Default Colormap 
    end
    
    % MatVTK configuration Parameters
    if ~isfield(config,'ambient')
        config.ambient            = 0;             % Ambient
    end
    if ~isfield(config,'ambientColor')
        config.ambientColor       = [1 1 1];       % Ambient Color
    end
    if ~isfield(config,'backgroundColor')
        config.backgroundColor    = [0 0 0];       % Background color
    end
    if ~isfield(config,'cameraFocalDisk')          % Camera Focal Disk
        config.cameraFocalDisk    = 1;
    end
    if ~isfield(config,'cameraFocalPoint')
        config.cameraFocalPoint   = [0 0 0];       % Camera Focal Point
    end
    if ~isfield(config,'cameraPosition')
        config.cameraPosition     = [29.0283 26.0049 164.4381]; % Camera Position
    end
    if ~isfield(config,'cameraViewUp')           
        config.cameraViewUp       = [0 0 1];    % Camera View Up
    end
    if ~isfield(config,'color')
        config.color              = [1 1 1];      %  Surface Color
    end 
    if ~isfield(config,'diffuse')                 
        config.diffuse            = 1;            % Use diffuse
    end
    if ~isfield(config,'diffuseColor')
        config.diffuseColor       = [1 1 1];      % Diffuse Color
    end
    if ~isfield(config,'edgeColor')               
        config.edgeColor          = [0 0 0];      % Edge Color
    end
    if ~isfield(config,'edgeVisibility')
        config.edgeVisibility     = 0;            % Edge Visibility
    end
    if ~isfield(config,'interpolation')
        config.interpolation      = 'gouroud';    % Interpolation
    end
    if ~isfield(config,'lightFollowCamera')
        config.lightFollowCamera  = 1;            % Light Following the camera
    end
    if ~isfield(config,'lineWidth')
        config.lineWidth          = 1;            % Line Width
    end
    if ~isfield(config,'opacity')
        config.opacity            = 1;            % Opacity
    end
    if ~isfield(config,'orientation')
        config.orientation        = [0 0 0];      % Orientation
    end
    if ~isfield(config,'origin')
        config.origin             = [0 0 0];      % Origin
    end
    if ~isfield(config,'position')
        config.position           = [0 0 0];      % Position
    end
    if ~isfield(config,'representation')
        config.representation     = 'surface';    % Object Representation
    end
    if ~isfield(config,'scale')
        config.scale              = [1 1 1];      % Scale
    end
    if ~isfield(config,'shading')
        config.shading            = 0;            %  Shading
    end 
    if ~isfield(config,'specular')                % Specular percent of lighting
        config.specular           = 0.6440;     
    end
    if ~isfield(config,'specularColor')
        config.specularColor      = [0.80 0.80 0.80];    % Specular Color
    end
    if ~isfield(config,'specularPower')           % Specular Power
        config.specularPower      = 100;      
    end
    if ~isfield(config,'twoSidedLighting')
        config.twoSidedLighting   = 1;            % Two Side Lighting
    end
end
%=========================================================================%
tr = config.tr;
sa = config.sa;
cl = config.cl;
%=========================Main Program====================================%
switch lower(sa)
    case 'y'
        sa = 1;
    case 'n'
        sa = 0;
end
[tr] = sort(tr);
wh = whos('Tracts');
if (strcmp(wh.class,'struct'))|(strcmp(wh.class,'cell'))
    Ns = size(Tracts,1);
elseif ischar(Tracts(1,:));
    Ns = size(Tracts,1);
end
tr((size(tr,2)+1):Ns) = max(tr);
tr = tr(1:Ns);


%col = [213 221 227;175 206 227;149 196 228; 120 186 230;87 172 231;24 146 232;6 73 119;244 207 154;244 192 117;244 179 86;244 161 43; 212 133 20;158 101 19; 113 71 12]/255;
col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
cont = 0;
Osys = computer;
for j = 1:Ns
    config.opacity = tr(j);
    if ischar(Tracts(1,:));
        [header tracks] = trk_read(deblank(Tracts(j,:)));
    end
    if strcmp(wh.class,'struct');
        tracks = Tracts;
    end
    if strcmp(wh.class,'cell');
        tracks = Tracts{j,1};
    end
    re = floor(Ns/Ncolor); col = repmat(col,[re+1 1]);
    if (sa ==1)&(j==1)
         VTK = vtkinit;
    elseif (sa ==0)
         VTK = vtkinit;
    end
    Points = [0 0 0];
    for  i = 1:size(tracks,2);
        Points = [Points;tracks(i).matrix];
        index(i) = size(Points,1);
    end
    Points(1,:) = [];
    index = [1 index];
    index(end) = [];
    config.color = col(j,:);
    config.diffuseColor = col(j,:);
    vtkplotlines(gcvtk, Points,index(:), config);
end
vtkorientation();
%========================End of main program==============================%
return

%=======================Internal functions================================%

function [trn,Surfa] = Reord_Tr(SurfF,tr,Ns,wh);

if strcmp(wh.class,'struct')
    Surfa{1,1} = SurfF;
elseif strcmp(wh.class,'cell')
    Surfa = SurfF;
end

for j = 1:Ns
    if ischar(SurfF(1,:));
        [ pth nm ext] = fileparts(SurfF(j,:));
        if strcmp(ext(1:4),'.mat');
            Surf = load('-mat',[pth filesep nm ext(1:4)]);
            Surfa{j,1} = Surf.Surf;
        else
            [OutFiles, Surf] = Exp_Surf(SurfF(j,:), '0', '','', 'imp','n');
            Surfa{j,1} = Surf{1};
        end
    end
    Surf = Surfa{j,1};
    for i = 1:size(Surf,2);
nam(j) =abs(min(Surf(i).SurfData.vertices(:,1)));
    end
end
[b,t] = sort(nam);
trn(t(end:-1:1)) = tr;

return
%=========================================================================%
%=========================Internal Functions===============================





function Plot_Surf_VTK(Surfa,config);
%
% Syntax :
% Plot_Surf_VTK(Surfa,tr,sa);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   Surfa       : Surfaces files.
%
%   config      : Configuration Structure
%         % Plot_Surf Options   
%         config.tr                 : Transparency vector(values have to be between 0 and 1).
%         config.sa                 : Boolean variable(config.sa = 1, plot surfaces in the same window,
%                                     config.sa = 0, plot surfaces in diferent windows)
%         config.cl                 : Colormap used to see the results.
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
%         config.triangulationTolerance = 1; % Triangulation Tolerance
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
    [Surfa,sts] = spm_select([1 Inf],'any','Selecting Surface Files','',cd);
end

if nargin <2
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
        config.cameraPosition     = [-160 0 0 ];  % Camera Position
        config.cameraViewUp       = [0 0 0];      % Camera View Up
        config.color              = [1 1 1];      % Surface Color
        config.diffuse            = 1;            % Diffuse
        config.diffuseColor       = [1 1 1];      % Diffuse Color
        config.edgeColor          = [0 0 0];      % Edge Color
        config.edgeVisibility     = 0;            % Edge Visibility
        config.interpolation      = 'phong';    % Interpolation
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
        config.triangulationTolerance = 0.1;        % Triangulation Tolerance
else
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
        config.cameraPosition     = [-160 0 0]; % Camera Position
    end
    if ~isfield(config,'cameraViewUp')           
        config.cameraViewUp       = [0 0 0];    % Camera View Up
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
        config.interpolation      = 'flat';       % Interpolation
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
    if ~isfield(config,'triangulationTolerance')
        config.triangulationTolerance = 0.1;        % Triangulation Tolerance
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
wh = whos('Surfa');
if (strcmp(wh.class,'struct'))|(strcmp(wh.class,'cell'))
    Ns = size(Surfa,1);
elseif ischar(Surfa(1,:));
    Ns = size(Surfa,1);
end
tr((size(tr,2)+1):Ns) = max(tr);
tr = tr(1:Ns);
[tr,Surfa] = Reord_Tr(Surfa,tr,Ns,wh);

%col = [213 221 227;175 206 227;149 196 228; 120 186 230;87 172 231;24 146 232;6 73 119;244 207 154;244 192 117;244 179 86;244 161 43; 212 133 20;158 101 19; 113 71 12]/255;
col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
cont = 0;
Osys = computer;
for j = 1:Ns
    config.opacity = tr(j);
    Surf = Surfa{j,1};
    re = floor(length(Surf)/Ncolor); col = repmat(col,[re+1 1]);
    if ~isfield(Surf,'SurfData');
        errordlg('This file does not contain surface information');
        continue;
    end
    if (sa ==1)&(j==1)
        if ~isfield(Surf(1),'Name')
            Surf(1).Name = 'Unknown';
        end
        nm = Surf(1).Name;
        VTK = vtkinit;
        
        %colordef black;h = figure('numbertitle','off','Color', 'black','name',['IBASPM Surface Ploting...:  Plot ' nm]);hold on;  grid on; box on;
    elseif (sa ==0)
        if ~isfield(Surf(1),'Name')
            Surf(1).Name = 'Unknown';
        end
        nm = Surf(1).Name;
        VTK = vtkinit;
        %colordef black;h = figure('numbertitle','off','Color','black','name',['IBASPM Surface Ploting...:  Plot ' nm]);hold on;  grid on; box on;
    end
    if size(Surf,2)~=1
        for i = 1:size(Surf,2);
            if isfield(Surf(i).SurfData,'VertexNormals');
                Surf(i).SurfData = rmfield(Surf(i).SurfData,'VertexNormals');
            end
            if isfield(Surf(i).SurfData,'faces')
                if isfield(Surf(i).SurfData,'FaceVertexCData')
                    Surf(i).SurfData.FaceColor = 'interp';
                    [a,b,c] = unique(Surf(i).Is);
                    Colors = [a Surf(i).SurfData.FaceVertexCData(b,:)];
                    strsurf = vtkplotmesh(gcvtk, Surf(i).SurfData.vertices, Surf(i).SurfData.faces, Surf.Is, Colors, config);
                else
                    strsurf = vtkplotmesh(gcvtk, Surf(i).SurfData.vertices, Surf(i).SurfData.faces, ones(max(Surf(i).SurfData.faces(:)),1), [1 col(i,:)], config);
                end
            end
        end
    else
        if isfield(Surf.SurfData,'VertexNormals');
            Surf.SurfData = rmfield(Surf.SurfData,'VertexNormals');
        end
        if ~isfield(Surf,'Is')&~isfield(Surf.SurfData,'FaceVertexCData')
            if isfield(Surf.SurfData,'faces')
                strsurf = vtkplotmesh(gcvtk, Surf.SurfData.vertices, Surf.SurfData.faces, ones(max(Surf.SurfData.faces(:)),1), [1 col(j,:)], config);
            end
        elseif isfield(Surf,'Is')&~isfield(Surf.SurfData,'FaceVertexCData')
            %% Plot using vertices
            [Colors] = Surf_Color(Surf,cl);
            Surf.SurfData.FaceVertexCData = Colors;
            Surf.SurfData.FaceColor = 'interp';
            [a,b,c] = unique(Surf.Is);
            Colors = [a Surf.SurfData.FaceVertexCData(b,:)];
            if isfield(Surf.SurfData,'faces')
                strsurf = vtkplotmesh(gcvtk, Surf.SurfData.vertices, Surf.SurfData.faces, Surf.Is, Colors, config);
            end
            if sum(Surf.Is-floor(Surf.Is)) ~=0
                h = colorbar;
                if ~strcmp(cl,'spectral')
                    colormap(cl);
                else
                    colormap(spectral(64));
                end
                range = max(Surf(1).Is)-min(Surf(1).Is);values =  min(Surf(1).Is):range/10: max(Surf(1).Is);
                set(h,'YTickLabel',num2str(values'));
            end
            Surf.SurfData = rmfield(Surf.SurfData,'FaceVertexCData');
        elseif ~isfield(Surf,'Is')&(isfield(Surf.SurfData,'FaceVertexCData'))
            Surf.SurfData.FaceColor = 'interp';
            [a,b,c] = unique(Surf.SurfData.FaceVertexCData,'rows');
            d = [1:length(a)]';
            Surf.Is = d(c);
            Colors = [d a];
            if isfield(Surf.SurfData,'faces')
                strsurf = vtkplotmesh(gcvtk, Surf.SurfData.vertices, Surf.SurfData.faces, Surf.Is, Colors, config);
            end
        end
    end
    if isfield(Surf(1),'Is')&isfield(Surf(1).SurfData,'FaceVertexCData')
        cont = cont + 1;
        Surf.SurfData.FaceColor = 'interp';
        [a,b,c] = unique(Surf.Is);
        Colors = [a Surf.SurfData.FaceVertexCData(b,:)];
        strsurf = vtkplotmesh(gcvtk, Surf.SurfData.vertices, Surf.SurfData.faces, Surf.Is, Colors, config);
        if sum(Surf(1).Is-floor(Surf(1).Is)) ~=0
            h = colorbar;
            if ~strcmp(cl,'spectral')
                colormap(cl);
            else
                colormap(spectral(64));
            end
            range = max(Surf(1).Is)-min(Surf(1).Is);values =  min(Surf(1).Is):range/10: max(Surf(1).Is); 
            set(h,'YTickLabel',num2str(values'));
        end
    end
    vtkorientation();
    vtkview(gcvtk,0,270)
    vtkview(gcvtk,-90,0)
end

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





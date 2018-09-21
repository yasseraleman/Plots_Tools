function hf = Plot_Graph_Connectivity_VTK(ConnectFile,config, plotopts);
%
% Syntax :
%  hf = Plot_Graph_Connectivity_VTK(ConnectFile);
%
% This script plots the Graph  %% IT NEEDS TO BE FINISHED
%
% Input Parameters:
%      ConnectFile      :  Connectivity Matrix
%       Graphprop       :  Graph Property (clustering, pathlength, etc)
% Output Parameters:
%          hf           :  Figure Handle
%
% Related references:
%
%
% See also:
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% May 18th 2012
% Version $1.0
close all;


%% =====================Checking Input Parameters======================= %%

%=====================Checking Input Parameters===========================%
ConnectFile = '/media/Data/PROCESSING_RESULTS/PEPS/7-connectome/0003-20090713/tractres/probtrack/0003-20090713-Connectivity_Matrix-aparc+aseg.txt';
% if nargin==0
%     [ConnectFile,sts] = spm_select([1 Inf],'any','Selecting Connectivity Files','',cd);
% end

if nargin <2
        % Point Options
        config.pointSize          = 3;             % The radius of the spheres can be set. (default=1.0)
        config.phiResolution      = 12;            % Number of latitude lines
        config.thetaResolution    = 12;            % Number of longitude lines

        % Lines Options
        config.lineWidth          = 1;             % Sets the line width, defaults to 1.0
        config.lineTubes          = 1;             % 0=draw lines, 1=draw tubes. defaults to 0.
        config.tubeSides          = 6;             % Number of sides for the tube, defaults to 6 (hexagonal surface) 
        config.tubeRadius         = 1;             % Tube thickness, defaults to 1.
    
        % Plot_Surf Options   
        config.tr                 = 1;            % Opacity vector
        config.sa                 = 'y';          % Plotting in the same window
        config.cl                 = 'jet';        % Default Colormap   
        
        % MatVTK Options
        config.ambient            = 0;            % Ambient
        config.ambientColor       = [1 1 1];      % Ambient Color
        config.backgroundColor    = [1 1 1];      % Background color
        config.cameraFocalDisk    = 1;            % Camera Focal Disk
        config.cameraFocalPoint   = [0 0 0];      % Camera Focal Point
        config.cameraPosition     = [0 0 160]; % Camera Position
        config.cameraViewUp       = [1 0 0];    % Camera View Up
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

     % Points Options
    if ~isfield(config,'pointSize')
        config.pointSize          = 3;  % The radius of the spheres can be set. (default=1.0)
    end
    if ~isfield(config,'phiResolution')
        config.phiResolution          = 12;  % Number of latitude lines
    end
    if ~isfield(config,'thetaResolution')
        config.thetaResolution          = 12;  % Number of longitude lines
    end
    
    % Lines Options
    if ~isfield(config,'lineWidth')
        config.lineWidth          = 1;  % Sets the line width, defaults to 1.0
    end
    if ~isfield(config,'lineTubes')
        config.lineTubes          = 1;  % 0=draw lines, 1=draw tubes. defaults to 0.
    end
    if ~isfield(config,'tubeSides')
        config.tubeSides          = 6;  % Number of sides for the tube, defaults to 6 (hexagonal surface)
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
        config.cameraPosition     = [0 0 160]; % Camera Position
    end
    if ~isfield(config,'cameraViewUp')           
        config.cameraViewUp       = [1 0 0];    % Camera View Up
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

if nargin <3
    plotopts.thr = 0.3;
    plotopts.bin = 0;
    plotopts.cl ='jet';
    plotopts.linethick = [1 config.pointSize-0.5]; % LIne Scaling
    plotopts.graphprop = 'clustering';
    plotopts.linesuniqcolor = 0; 
end
if ~isfield('plotopts','thr')
    plotopts.thr = 0.3;
end
if ~isfield('plotopts','bin')
    plotopts.bin = 0;
end
if ~isfield('plotopts','cl')
    plotopts.cl ='jet';
end
if ~isfield('plotopts','scalLine')
    plotopts.scalLine = 8; % LIne Scaling
end
if ~isfield('plotopts','scalPlot')
    plotopts.scalPlot = 12; % LIne Scaling
end
if ~isfield('plotopts','graphprop')
    plotopts.graphprop = 'clustering';
end
if ~isfield('plotopts','atlastype')
    plotopts.atype = 'aparc+aseg';
end
if ~isfield('plotopts','linesuniqcolor')
    plotopts.linesuniqcolor = 1;
end
if ~isfield('plotopts','pointsuniqcolor')
    plotopts.pointsuniqcolor = 0;
end


%% ===============End of Checking Input Parameters====================== %%

%% ================= Reading the Connectivity Matrix =================== %%
if ischar(ConnectFile)
    Connect = Read_Connectivity_Matrix(ConnectFile,0);
    Connect.Matrix = Connect.Matrix/max(Connect.Matrix(:));
else
    Connect.Matrix = ConnectFile;
    if ~isfield('plotopts','atlastype')
        plotopts.atype = 'aparc+aseg';
    end
end

%% ============ End of Reading the Connectivity Matrix ================= %%

%% ============ Thresholding the Connectivity Matrix =================== %%
ind = find(Connect.Matrix<=plotopts.thr);
Connect.Matrix(ind) = 0;
if plotopts.bin
    Connect.Matrix =  logical(Connect.Matrix);
end
Sizes = ones(size(Connect.Matrix,1),1)*12;
%% ============ End of Thresholding the Connectivity Matrix ============ %%
% ------------------- Detecting the Atlas Type -------------------------- %
if ischar(ConnectFile)
    fid  = fopen(ConnectFile,'rt');
    linet = fgetl(fid);linet = fgetl(fid);linet = fgetl(fid);
    [temp,atype] = strread(linet,'%s%s','delimiter',':');
    atype = char(atype);
    fclose(fid);
else
    atype = plotopts.atype;
end
[SIds,Names,Colors] = Gray_Matter_codes(atype);
if size(Connect.Matrix,1) == 68
    Names = Names(15:end,:);
    Colors = Colors(15:end,:);
end
Colors = Colors./255;
if nargin > 2
    switch plotopts.graphprop
        case 'clustering'
            txt = clustering_coef_wu(Connect.Matrix);
        case 'pathlength'
            
        otherwise
    end
    Colors = Val2colors(txt,plotopts.cl);
    mini = min(nonzeros(txt(:)));
    maxi = max(nonzeros(txt(:)));
    Sizes = (txt - mini).*(plotopts.scalPlot - 3)./(maxi - mini) + 3;
end

% ------------------- Loading The Mass Center Coordinates --------------- %
CoordFile = which([atype '-MassCenterCoords.txt']);
%[temp,Xc,Yc,Zc] = textread(CoordFile,'%s%f%f%f','delimiter',';');

[Xc,Yc,Zc] = textread('/media/COSAS/scripts/Morpho_Stuff/Centers.txt','%f%f%f','delimiter',' ');
if size(Connect.Matrix,1) == 68
    Xc = Xc(15:end);
    Yc = Yc(15:end);
    Zc = Zc(15:end);
end
%% ============ Selecting The Upper triangle in Connectivity Matrix ==== %%
Tmask = triu(Connect.Matrix);
%% ====== End of Selecting The Upper triangle in Connectivity Matrix ==== %
load('/media/COSAS/scripts/Morpho_Stuff/SurfacesJ.mat');
config.cameraFocalPoint =  [mean(Xc) mean(Yc) mean(Zc)];
%==========================  Main Program   ==============================%

%% ===================== Plotting Connections =========================== %
% --------- Line Colors according to the orientation of connections ----- %
LabNet = SubNetworks_Extraction(Connect.Matrix,0);
ind = find(Tmask);
[ix,iy] = ind2sub(size(Tmask),ind);
NetId = unique(LabNet(ind));

% -------------------- Lines Colors ---------------------------------- %
col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
re = floor(length(NetId)/Ncolor); col = repmat(col,[re+1 1]);
ColorsL = col(LabNet(ind),:);

% ---------- Line Width according to Connectivity values ---------------- %
[ut,i1t,it2] = unique(Connect.Matrix(ind));
LG = linspace(plotopts.linethick(1),plotopts.linethick(2),length(ut))';
LG = LG(it2);

vtkinit;
% ------------------------------ Plotting Surfaces ---------------------- %
if config.backgroundColor == [0 0 0];
    config.color = [1 1 1];
    if plotopts.linesuniqcolor
        ColorsL = repmat([ 0    0.3804    0.5961],[size(ColorsL,1) 1]);
    end
    if plotopts.pointsuniqcolor
        Colors = repmat([0.9961    0.9451    0.3804],[size(Colors,1) 1]);
    end
elseif config.backgroundColor == [1 1 1];
    config.color = [0 0 0];
    if plotopts.linesuniqcolor
        ColorsL = repmat([0.9961    0.9451    0.3804],[size(ColorsL,1) 1]);
    end
    if plotopts.pointsuniqcolor
        Colors = repmat([ 0    0.3804    0.5961],[size(Colors,1) 1]);
    end
end

% ------------------------------ Plotting Mesh ------------------------ %
temp = config.opacity;
config.opacity = 0.03;
config.origin  = [mean(Xc) mean(Yc) mean(Zc)];
strsurf = vtkplotmesh(gcvtk, Surfn.SurfData.vertices, Surfn.SurfData.faces, config);
config.opacity = temp;


% ------------------------------ Plotting Points ------------------------ %
indexes = unique([ix;iy]);
vtkplotpoints([Xc(indexes) Yc(indexes) Zc(indexes)], [1:length(indexes)]', [[1:length(indexes)]' Colors(indexes,:)], config);
indexesr = find(ismember([1:length(Xc)],indexes) == 0);
temp = config.opacity;
config.opacity = 0.1;
config.pointSize = 1.5;
vtkplotpoints([Xc(indexesr) Yc(indexesr) Zc(indexesr)], [1:length(indexesr)]', [[1:length(indexesr)]' repmat([0.4 0.4 0.4],[length(indexesr) 1])], config);
config.opacity = temp;
% ------------------------------ Plotting Lines ------------------------ %
for i = 1:length(ind)
    config.color = ColorsL(i,:);
    config.diffuseColor = ColorsL(i,:);
    config.lineWidth = LG(i);
    vtkplotlines(gcvtk, [Xc(ix(i)) Yc(ix(i)) Zc(ix(i));Xc(iy(i)) Yc(iy(i)) Zc(iy(i))], config);
end
vtkorientation(gcvtk,config);
%% ===================== End of Plotting Connections ==================== %
return;

function Colors = Val2colors(txt,cl);
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

function Names = Load_Struct_Names;

Names = strvcat('Left-Thalamus-Proper',...
    'Left-Caudate',...
    'Left-Putamen',...
    'Left-Pallidum',...
    'Left-Hippocampus',...
    'Left-Amygdala',...
    'Left-Accumbens-area',...
    'Right-Thalamus-Proper',...
    'Right-Caudate',...
    'Right-Putamen',...
    'Right-Pallidum',...
    'Right-Hippocampus',...
    'Right-Amygdala',...
    'Right-Accumbens-area',...
    'ctx-lh-bankssts',...
    'ctx-lh-caudalanteriorcingulate',...
    'ctx-lh-caudalmiddlefrontal',...
    'ctx-lh-cuneus',...
    'ctx-lh-entorhinal',...
    'ctx-lh-fusiform',...
    'ctx-lh-inferiorparietal',...
    'ctx-lh-inferiortemporal',...
    'ctx-lh-isthmuscingulate',...
    'ctx-lh-lateraloccipital',...
    'ctx-lh-lateralorbitofrontal',...
    'ctx-lh-lingual',...
    'ctx-lh-medialorbitofrontal',...
    'ctx-lh-middletemporal',...
    'ctx-lh-parahippocampal',...
    'ctx-lh-paracentral',...
    'ctx-lh-parsopercularis',...
    'ctx-lh-parsorbitalis',...
    'ctx-lh-parstriangularis',...
    'ctx-lh-pericalcarine',...
    'ctx-lh-postcentral',...
    'ctx-lh-posteriorcingulate',...
    'ctx-lh-precentral',...
    'ctx-lh-precuneus',...
    'ctx-lh-rostralanteriorcingulate',...
    'ctx-lh-rostralmiddlefrontal',...
    'ctx-lh-superiorfrontal',...
    'ctx-lh-superiorparietal',...
    'ctx-lh-superiortemporal',...
    'ctx-lh-supramarginal',...
    'ctx-lh-frontalpole',...
    'ctx-lh-temporalpole',...
    'ctx-lh-transversetemporal',...
    'ctx-lh-insula',...
    'ctx-rh-bankssts',...
    'ctx-rh-caudalanteriorcingulate',...
    'ctx-rh-caudalmiddlefrontal',...
    'ctx-rh-cuneus',...
    'ctx-rh-entorhinal',...
    'ctx-rh-fusiform',...
    'ctx-rh-inferiorparietal',...
    'ctx-rh-inferiortemporal',...
    'ctx-rh-isthmuscingulate',...
    'ctx-rh-lateraloccipital',...
    'ctx-rh-lateralorbitofrontal',...
    'ctx-rh-lingual',...
    'ctx-rh-medialorbitofrontal',...
    'ctx-rh-middletemporal',...
    'ctx-rh-parahippocampal',...
    'ctx-rh-paracentral',...
    'ctx-rh-parsopercularis',...
    'ctx-rh-parsorbitalis',...
    'ctx-rh-parstriangularis',...
    'ctx-rh-pericalcarine',...
    'ctx-rh-postcentral',...
    'ctx-rh-posteriorcingulate',...
    'ctx-rh-precentral',...
    'ctx-rh-precuneus',...
    'ctx-rh-rostralanteriorcingulate',...
    'ctx-rh-rostralmiddlefrontal',...
    'ctx-rh-superiorfrontal',...
    'ctx-rh-superiorparietal',...
    'ctx-rh-superiortemporal',...
    'ctx-rh-supramarginal',...
    'ctx-rh-frontalpole',...
    'ctx-rh-temporalpole',...
    'ctx-rh-transversetemporal',...
    'ctx-rh-insula');
return
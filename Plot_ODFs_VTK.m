function OutFigures = Plot_ODFs_VTK(BG_Image, ODF_Image, ODF_Dirs, slices, plotopts);
%
% Syntax :
% Images = Image_plus_overlay(BG_Image,OverlayImage,slices);
%
% The following program shows how to view brain slices and to overlay maps
% on top of anatomical scans.
%
% Input Parameters:
%       BG_Image       : BackGround Image
%       InputImage     : ODF Image
%       ODF_Dirs       : ODF Directions
%        slices        : Nx3 vector specifying wich slices will be ploted.
%                         ie. [Sagital_Slice Coronal_Slice  Axial_Slice];
%        plotopts      :  Plot options
%                          plotopts.transp         : Overlay Transparency (0 to 1)
%                          plotopts.lims.Xlim      : [Xmin Xmax] X Limits (voxels)
%                          plotopts.lims.Ylim      : [Ymin Ymax] Y Limits (voxels)
%                          plotopts.lims.Zlim      : [Zmin Zmax] Z Limits (voxels)
%
%  Output Parameters:
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
%
BG_Image = '/media/COSAS/Test/Erick/gfa_small.nii';
%  % BG_Image = '/media/Data/PROCESSING_RESULTS/ASPERGER/5-freesurfer_processing/ASPER_00001__101-20060510/tmp/T1t.nii'
ODF_Image = '/media/COSAS/Test/Erick/ODF_matrix_4D_small.nii';
% % % %OverlayImage = '/home/yaleman/Data/5-freesurfer_processing/ASPER_00001__101-20060510/tmp/aparc+asegt.nii'
%  OverlayImage = '/media/Data/PROCESSING_RESULTS/ASPERGER/7-connectome/ASPER_00001__101-20060510/dtreconst/dtk/ASPER_00001__101-20060510_fa.nii'
% % % %OverlayImage = '/home/yaleman/Data/7-connectome/ASPER_00001__101-20060510/dtreconst/fsl/ASPER_00001__101-20060510_fa_color.nii'
% % % OverlayImage = '/home/yaleman/Data/7-connectome/ASPER_00001__101-20060510/dtreconst/dtk/ASPER_00001__101-20060510_v1.nii'
slices = [0 1 0];

ODF_Dirs = '/media/COSAS/Test/Erick/724_shell.txt';


%% ================== Cheking Input Parameters ========================= %%
close all;
V = spm_vol(BG_Image);
VO = spm_vol(deblank(ODF_Image(1,:)));
[pth,name,ext] = fileparts(V(1).fname);
[ptho,nameo,exto] = fileparts(VO(1).fname);

if nargin < 4
    plotopts.colormap = 'jet';
    %          plotopts.lims.Xlim = [30 70];
    %          plotopts.lims.Ylim = [40 70];
    %     plotopts.lims.Zlim = [10 50];
end

dirs = load(deblank(ODF_Dirs));
faces = convhulln(dirs);
% dirs = [dirs;-1*dirs];

if ~isfield(plotopts,'lims')
    plotopts.lims.Xlim = [1 V(1).dim(1)];
    plotopts.lims.Ylim = [1 V(1).dim(2)];
    plotopts.lims.Zlim = [1 V(1).dim(3)];
end
if ~isfield(plotopts.lims,'Xlim')
    plotopts.lims.Xlim = [1 V(1).dim(1)];
end

if ~isfield(plotopts.lims,'Ylim')
    plotopts.lims.Ylim = [1 V(1).dim(2)];
end

if ~isfield(plotopts.lims,'Zlim')
    plotopts.lims.Zlim = [1 V(1).dim(3)];
end

%% ============== End of Cheking Input Parameters ====================== %%
voxsize = sqrt(sum(V(1).mat(1:3,1:3).^2));
% TitleF = ['Image + Overlay: ' upper(name) ' + ' upper(nameo)];
% ATtile = [ upper(name) '_plus_' upper(nameo)];
% OutFigures = '';
for t = 1:size(slices,1)
    slice = slices(t,:);
    ind = find(slice);
    for i = 1:length(ind)
%         colordef black;
%         hf = figure('numbertitle','off','name',TitleF,'Color',[0 0 0],'Position',[0 0 1200 900]);
        
        %% ============ Plotting ODFs ========================== %%
        if ind(i) == 3
            DIM = V.dim([1 2]);
            C=[1 0 0 0;0 1 0 0;0 0 1 slice(3);0 0 0 1]; % Axial Transformation Matrix
            Xlim = [plotopts.lims.Xlim(1):plotopts.lims.Xlim(2)]; % Axial Limits
%             Ylim = [plotopts.lims.Ylim(1):plotopts.lims.Ylim(2)];
%             nametotal = [ATtile '-Axial_' num2str(slice(3))];
%             DTtitle = [ATtile ' : Axial Slice Number ' num2str(slice(3)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Y from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
        elseif ind(i) == 2
            DIM = V.dim([1 3]);
            C=[1 0 0 0;0 0 1 slice(2);0 1 0 0;0 0 0 1]; % Coronal Transformation Matrix
            Xlim = [plotopts.lims.Xlim(1):plotopts.lims.Xlim(2)]; % Coronal Limits
            Ylim = [plotopts.lims.Zlim(1):plotopts.lims.Zlim(2)];
%             nametotal = [ATtile '-Coronal_' num2str(slice(2))];
%             DTtitle = [ATtile ' : Coronal Slice Number ' num2str(slice(2)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
        elseif ind(i) == 1
            DIM = V.dim([2 3]);
            C=[0 0 1 slice(1);1 0 0 0;0 1 0 0;0 0 0 1]; % Sagital Transformation Matrix
            Xlim = [plotopts.lims.Ylim(1):plotopts.lims.Ylim(2)]; % Sagital Limits
            Ylim = [plotopts.lims.Zlim(1):plotopts.lims.Zlim(2)];
%             nametotal = [ATtile '_Sag_' num2str(slice(1))];
%             DTtitle = [ATtile ' : Sagital Slice Number ' num2str(slice(1)) '. Limits: Y from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
        end
        
        % Reading Background Image
        Is = spm_slice_vol(V,C,DIM,0);
        Is = Is(Xlim,Ylim);
%         Is=rot90(Is);
%         imagesc(Is);
%         colormap gray;axis equal;axis off;
%         hold on;
        for k = 1:length(VO)
            IOst = spm_slice_vol(VO(k),C,DIM,0);
            indd = isnan(IOst);
            IOst(indd) = 0;
            IOst = IOst(Xlim,Ylim);
            IOs(:,:,k) = IOst;
%             IOs(:,:,k) = rot90(IOst);
        end
        indd = find((sum(IOs,3)~=0)&(sum(IOs,3)~=NaN));
        [X,Y] = ind2sub([size(IOs,1) size(IOs,2)],indd);
        % ============ Creating Patch ============================%
        verts = [0 0 0]; F = [0 0 0]; Vt = [0 0 0];
        %vtkinit; % FOR VTK
        for k = 1:30%length(indd)
            k
            ODF = IOs(X(k),Y(k),:);ODF = [ODF(:);ODF(:)];
            Vs = (repmat((1/2)*ODF./max(ODF),[1 3]).*dirs);
            Vt = [Vt;Vs];
            verts = [verts;Vs+ repmat( [Y(k) X(k) 2],[size(dirs,1) 1])];
            F = [F;faces+max(F(:))];
           % a = Vs ./(repmat(sqrt(dot(Vs,Vs,2)),1,3) + eps);
           % strsurf = vtkplotmesh(gcvtk, Vs+ repmat( [Y(k) X(k) 0],[size(dirs,1) 1]), faces, [1:max(faces(:))]', [[1:max(faces(:))]' abs(a)]);
        end
        Vt(1,:) = [];verts(1,:) = []; F(1,:) = [];
        Surf.SurfData.vertices = verts;
        Surf.SurfData.faces = F;
        a = Vt ./(repmat(sqrt(dot(Vt,Vt,2)),1,3) + eps);
        Surf.SurfData.FaceVertexCData =abs(a);
        %vtktrisurf(Surf.SurfData.faces,Surf.SurfData.vertices(:,1),Surf.SurfData.vertices(:,2),Surf.SurfData.vertices(:,3),'ODFs',[pwd filesep 'vtktrisurf.vtu']);
        disp(['%%%%% File Save To ' pwd filesep 'vtktrisurf.vtu. Open ParaView to see the results']);
         Plot_Surf_VTK(Surf); % Plot with VTK
        
        % Plot with Matlab
%         patch('Vertices', Surf.SurfData.vertices, 'Faces', Surf.SurfData.faces,'FaceVertexCData', Surf.SurfData.FaceVertexCData,'FaceColor', 'interp', 'FaceLighting','phong',...
%             'LineStyle','-',...
%             'LineWidth',.05,...
%             'EdgeColor',[.3 .3 .3],...
%             'AmbientStrength',.4,...
%             'FaceLighting','phong',...
%             'SpecularColorReflectance',.2,...
%             'DiffuseStrength',.5,...
%             'BackFaceLighting', 'reverselit');
%         camlight;
%         axis equal
        %% ========== End of Plotting ODFs ======================%%
    end
    clear IOs;
    %% ============== End of Working With 4D Maps=============== %%
    
end
% ------------------ Removing Temporary Files ----------------------------%
% delete(deblank(ODF_Image(1,:)));
% delete(deblank(ODF_Image(1,:)));
% ------------------End of  Removing Temporary Files ---------------------%
return;

function Filenames = change_spacen_custom(Ii0,If0,order,Filenames);
if size(If0,1) == 1,
    If0 = If0(ones(size(Ii0,1),1),:);
elseif ~(size(Ii0,1)==size(If0,1)),
    error('it must be the same number of images')
end
outnames = '';
for image = 1:size(Ii0),
    Ii = deblank(Ii0(image,:));
    If = deblank(If0(image,:));
    Vi = spm_vol(Ii); n = length(Vi); name = Vi.fname;
    Vf = spm_vol(If);
    [a1,a2,a3] = fileparts(Ii); [b1,b2,b3] = fileparts(If);
    counter = 0; %handle = waitbar(counter,'resampling...');
    [X,Y] = ndgrid(1:Vf.dim(1),1:Vf.dim(2));
    for i = 1:n,
        Vn = Vf; c = spm_bsplinc(Vi(i),[order*ones(1,3) 0 0 0]);
        Vn.pinfo = Vi(i).pinfo; Vn.dt = Vi(i).dt; Vn.n = Vi(i).n;
        [pathstr,name,ext] = fileparts(Vi(i).fname);
        if nargin<4
            Vn.fname = fullfile(pathstr,['s' name ext]);
        elseif nargin == 4
            Vn.fname = deblank(Filenames(image,:));
        end
        try, Vn.userdata.g = Vf.userdata.g; Vn.userdata.mat = Vn.mat; end
        Vn = spm_create_vol(Vn); M = Vi(i).mat\Vn.mat;
        Nslice = Vn.dim(3);
        for z = 1:Nslice,
            counter = counter+1;
            A = spm_bsplins(c,M(1,1)*X + M(1,2)*Y + M(1,3)*z + M(1,4),...
                M(2,1)*X + M(2,2)*Y + M(2,3)*z + M(2,4),...
                M(3,1)*X + M(3,2)*Y + M(3,3)*z + M(3,4),...
                [order*ones(1,3) 0 0 0]);
            spm_write_plane(Vn,A,z);
            disp(['Resampling volume ' num2str(i) ' of ' num2str(n) ': Slice ' num2str(z) ' of ' num2str(Nslice)]);
            %waitbar(counter/(Vn.dim(3)*n),handle)
        end
        outnames = strvcat(outnames,Vn.fname);
    end, %close(handle);
end
outnames = unique(outnames,'rows');
fclose all;
return

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
    config.specular           = 0;       % Specular percent of lighting
    config.specularColor      = [0 0 0];    % Specular Color
    config.specularPower      = 0;          % Specular power
    config.twoSidedLighting   = 1;            % Two Side Lighting
    %config.triangulationTolerance = 0.1;        % Triangulation Tolerance
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
%     if ~isfield(config,'triangulationTolerance')
%         config.triangulationTolerance = 0.1;        % Triangulation Tolerance
%     end
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





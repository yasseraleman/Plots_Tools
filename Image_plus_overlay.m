function OutFigures = Image_plus_overlay(InputImage,OverlayImage,slices, plotopts, OutDir);
%
% Syntax :
% Images = Image_plus_overlay(InputImage,OverlayImage,slices);
%
% The following program shows how to view brain slices and to overlay maps
% on top of anatomical scans.
%
% Input Parameters:
%     InputImage       : Input Anatomical Image
%       OverlayImage   : Overlay Image
%        slices        : Nx3 vector specifying wich slices will be ploted.
%                         ie. [Sagital_Slice Coronal_Slice  Axial_Slice];
%        plotopts      :  Plot options
%                          plotopts.atype          : Atlas Type ('aparc+aseg','a2009s+aseg',
%                                                                'a2005s+aseg','ibaspm116',
%                                                                'ibaspm71')
%                          plotopts.transp         : Overlay Transparency (0 to 1)
%                          plotopts.colormap       : Overlay Color map ('hot','hsv','pink', etc)
%                          plotopts.itype          : Image Type ('atlas', 'mask' or 'image')
%                          plotopts.dtype          : Dimension Type ('scalar' or 'vectorial')
%                          plotopts.dtype.modul    : Vector field modulation (1: Modulate, 0: Otherwise)
%                          plotopts.linelength     : Vector line length
%                          plotopts.lims.Xlim      : [Xmin Xmax] X Limits (voxels)
%                          plotopts.lims.Ylim      : [Ymin Ymax] Y Limits (voxels)
%                          plotopts.lims.Zlim      : [Zmin Zmax] Z Limits (voxels)
%         OutDir        : Output Directory
%
%  Output Parameters:
%       OutFigures      : Figures Filenames
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
%  InputImage = '/media/Data/PROCESSING_RESULTS/ASPERGER/7-connectome/ASPER_00001__101-20060510/dtreconst/dtk/ASPER_00001__101-20060510_b0.nii';
%  % InputImage = '/media/Data/PROCESSING_RESULTS/ASPERGER/5-freesurfer_processing/ASPER_00001__101-20060510/tmp/T1t.nii'
% % %  %OverlayImage = '/home/yaleman/Data/5-freesurfer_processing/ASPER_00001__101-20060510/tmp/brainmaskt.nii'
% % % %OverlayImage = '/home/yaleman/Data/5-freesurfer_processing/ASPER_00001__101-20060510/tmp/aparc+asegt.nii'
%  OverlayImage = '/media/Data/PROCESSING_RESULTS/ASPERGER/7-connectome/ASPER_00001__101-20060510/dtreconst/dtk/ASPER_00001__101-20060510_fa.nii'
% % % %OverlayImage = '/home/yaleman/Data/7-connectome/ASPER_00001__101-20060510/dtreconst/fsl/ASPER_00001__101-20060510_fa_color.nii'
% % % OverlayImage = '/home/yaleman/Data/7-connectome/ASPER_00001__101-20060510/dtreconst/dtk/ASPER_00001__101-20060510_v1.nii'
% slices = [64 64 30];

% OverlayImage = strvcat('/home/yaleman/Data/7-connectome/ASPER_00001__101-20060510/dtreconst/dtk/ASPER_00001__101-20060510_odf.nii',...
% '/home/yaleman/Data/7-connectome/ASPER_00001__101-20060510/dtreconst/dtk/ASPER_00001__101-20060510_odfdirs.txt',...
% '/home/yaleman/Data/7-connectome/ASPER_00001__101-20060510/dtreconst/dtk/ASPER_00001__101-20060510_odffaces.txt');


%% ================== Cheking Input Parameters ========================= %%
close all;
V = spm_vol(InputImage);
VO = spm_vol(deblank(OverlayImage(1,:)));
[pth,name,ext] = fileparts(V(1).fname);
[ptho,nameo,exto] = fileparts(VO(1).fname);
if nargin < 4
    plotopts.transp = 0.4;
    plotopts.colormap = 'jet';
    plotopts.itype = 'image';
    plotopts.atype = 'aparc+aseg';
    plotopts.dtype.type = 'vectorial';
    plotopts.dtype.modul = 0;
    plotopts.linelength = 1/2;
    plotopts.figurevisible = 'off';
    plotopts.figcolor = 'black';
    %     plotopts.lims.Xlim = [20 120];
    %     plotopts.lims.Ylim = [20 120];
    %     plotopts.lims.Zlim = [10 50];
    

    OutDir = pth;
end
if nargin < 5
    OutDir = pth;
end

if size(OverlayImage,1)>1
    dirs = load(deblank(OverlayImage(2,:)));
    dirs(2) = -1*dirs(2); % 180 degrees rotation around X axis
    faces = load(deblank(OverlayImage(3,:)));
end
if ~isfield(plotopts,'transp')
    plotopts.transp = .4;
end
if ~isfield(plotopts,'colormap')
    plotopts.colormap = 'jet';
end
if ~isfield(plotopts,'itype')
    plotopts.itype = 'image';
end
if ~isfield(plotopts,'atype')
    plotopts.atype = 'aparc+aseg';
end

if isfield(plotopts,'dtype')
    if ~isfield(plotopts.dtype,'type')
        plotopts.dtype.type = 'vectorial';
    end
    if ~isfield(plotopts.dtype,'modul')
        plotopts.dtype.modul = 0;
    end
    
end
if ~isfield(plotopts,'linelength')
    plotopts.linelength = .5;
end


if ~isfield(plotopts,'lims')
    plotopts.lims.Xlim = [1 V.dim(1)];
    plotopts.lims.Ylim = [1 V.dim(2)];
    plotopts.lims.Zlim = [1 V.dim(3)];
end
if ~isfield(plotopts.lims,'Xlim')
    plotopts.lims.Xlim = [1 V.dim(1)];
end

if ~isfield(plotopts.lims,'Ylim')
    plotopts.lims.Ylim = [1 V.dim(2)];
end
if ~isfield(plotopts.lims,'Zlim')
    plotopts.lims.Zlim = [1 V.dim(3)];
end
if ~isfield(plotopts,'figurevisible')
    plotopts.figurevisible = 'off';
end
if ~isfield(plotopts,'figcolor')
    plotopts.figcolor = 'black';
end



%% ============== End of Cheking Input Parameters ====================== %%
voxsize = sqrt(sum(V.mat(1:3,1:3).^2));
TitleF = ['Image + Overlay: ' upper(name) ' + ' upper(nameo)];
ATtile = [ upper(name) '_plus_' upper(nameo)];
OutFigures = '';
for t = 1:size(slices,1)
    slice = slices(t,:);
    ind = find(slice);
    for i = 1:length(ind)
        switch plotopts.figcolor
            case 'black'
                colordef black;
                colordef black; hf = figure('numbertitle','off','name',TitleF,'Color',[0 0 0],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
            case 'white'
                colordef white;
                colordef white; hf = figure('numbertitle','off','name',TitleF,'Color',[1 1 1],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
        end
        if length(VO) == 1 % Verifying if the image is Coloured or not
            %% Non-Coloured Images
            if ind(i) == 3
                DIM = V.dim([1 2]);
                C=[1 0 0 0;0 1 0 0;0 0 1 slice(3);0 0 0 1]; % Axial Transformation Matrix
                Xlim = [plotopts.lims.Xlim(1):plotopts.lims.Xlim(2)]; % Axial Limits
                Ylim = [plotopts.lims.Ylim(1):plotopts.lims.Ylim(2)];
                nametotal = [ATtile '-Axial_' num2str(slice(3))];
                DTtitle = [ATtile ' : Axial Slice Number ' num2str(slice(3)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Y from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
            elseif ind(i) == 2
                DIM = V.dim([1 3]);
                C=[1 0 0 0;0 0 1 slice(2);0 1 0 0;0 0 0 1]; % Coronal Transformation Matrix
                Xlim = [plotopts.lims.Xlim(1):plotopts.lims.Xlim(2)]; % Coronal Limits
                Ylim = [plotopts.lims.Zlim(1):plotopts.lims.Zlim(2)];
                nametotal = [ATtile '-Coronal_' num2str(slice(2))];
                DTtitle = [ATtile ' : Coronal Slice Number ' num2str(slice(2)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
            elseif ind(i) == 1
                DIM = V.dim([2 3]);
                C=[0 0 1 slice(1);1 0 0 0;0 1 0 0;0 0 0 1]; % Sagital Transformation Matrix
                Xlim = [plotopts.lims.Ylim(1):plotopts.lims.Ylim(2)]; % Sagital Limits
                Ylim = [plotopts.lims.Zlim(1):plotopts.lims.Zlim(2)];
                nametotal = [ATtile '-Sagital_' num2str(slice(1))];
                DTtitle = [ATtile ' : Sagital Slice Number ' num2str(slice(1)) '. Limits: Y from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
            end
            % Reading Background Image
            Is = spm_slice_vol(V,C,DIM,0);
            Is = Is(Xlim,Ylim);% Selecting Limits
            Is=rot90(Is);
            
            % Reading Overlay Image
            IOs = spm_slice_vol(VO,C,DIM,0);
            IOs = IOs(Xlim,Ylim);% Selecting Limits
            IOs=rot90(IOs);
            
            % Plot Background
            imagesc(Is);colormap gray;axis equal;axis off;
            %% ============== Is it a Mask? ================================ %%
            sts = unique(IOs(:));
            sts(sts==0) = [];
            %% ============== End of Is it a Mask? ========================= %%
            if length(sts)==1
                %% ================== Plot Masks =========================== %%
                hold on;
                red = cat(3, ones(size(IOs,1),size(IOs,2)), zeros(size(IOs,1),size(IOs,2)), zeros(size(IOs,1),size(IOs,2)));
                h = imshow(red);
                hold off;
                set(h, 'AlphaData', plotopts.transp*IOs);
                %% ============== End of Plot Masks ======================== %%
            else
                %% ============== Plot Color Maps ========================== %%
                switch plotopts.itype
                    case 'image'
                        Colors = Image2Colors(IOs,plotopts.colormap);
                        hold on;
                        h = imshow(Colors);
                        hold off;
                        indd = find(IOs);
                        IOst = IOs*0;
                        IOst(indd) = 1;
                        set(h, 'AlphaData', plotopts.transp*IOst);
                        clear IOs;
                        %% ============== End of Plot Color Maps =============== %%
                    case 'mask'
                        hold on;
                        red = cat(3, ones(size(IOs,1),size(IOs,2)), zeros(size(IOs,1),size(IOs,2)), zeros(size(IOs,1),size(IOs,2)));
                        h = imshow(red);
                        hold off;
                        indd = find(IOs);
                        IOst = IOs*0;
                        IOst(indd) = 1;
                        set(h, 'AlphaData', plotopts.transp*IOst);
                    otherwise
                        %% ================== Plot Atlases ===================== %%
                        [GMcodes,Names, Colors] = Brain_GM_codes(plotopts.atype);
                        sts = unique(IOs(:));
                        sts(sts==0) = [];
                        if isempty(GMcodes)
                            col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
                            Ncolor = size(col,1);
                            re = floor(length(sts)/Ncolor); col = repmat(col,[re+1 1]);
                            Colors = col(1:length(sts),:)*255;
                            GMcodes = sts;
                        end
                        indd = find(ismember(IOs(:),GMcodes) == 0);
                        IOs(indd) = 0;
                        sts = unique(IOs(:));
                        sts(sts==0) = [];
                        for j = 1:length(sts)
                            IOst = IOs*0;
                            indc = find(GMcodes == sts(j));
                            inds = find(IOs == sts(j));
                            IOst(inds) = 1;
                            color = Colors(indc,:)/255;
                            hold on;
                            red = cat(3, ones(size(IOs,1),size(IOs,2))*color(1), ones(size(IOs,1),size(IOs,2))*color(2), ones(size(IOs,1),size(IOs,2))*color(3));
                            h = imshow(red);
                            hold off;
                            set(h, 'AlphaData', plotopts.transp*IOst);
                        end
                        %% =============== End of Plot Atlases ================= %%
                end
            end
        else
            %% ============== Working With 4D Maps ========================= %%
            if  length(VO) == 3
                if ind(i) == 3
                    DIM = V.dim([1 2]);
                    C=[1 0 0 0;0 1 0 0;0 0 1 slice(3);0 0 0 1]; % Axial Transformation Matrix
                    Xlim = [plotopts.lims.Xlim(1):plotopts.lims.Xlim(2)]; % Axial Limits
                    Ylim = [plotopts.lims.Ylim(1):plotopts.lims.Ylim(2)];
                    nametotal = [ATtile '-Axial_' num2str(slice(3))];
                    DTtitle = [ATtile ' : Axial Slice Number ' num2str(slice(3)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Y from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
                elseif ind(i) == 2
                    DIM = V.dim([1 3]);
                    C=[1 0 0 0;0 0 1 slice(2);0 1 0 0;0 0 0 1]; % Coronal Transformation Matrix
                    Xlim = [plotopts.lims.Xlim(1):plotopts.lims.Xlim(2)]; % Coronal Limits
                    Ylim = [plotopts.lims.Zlim(1):plotopts.lims.Zlim(2)];
                    nametotal = [ATtile '-Coronal_' num2str(slice(2))];
                    DTtitle = [ATtile ' : Coronal Slice Number ' num2str(slice(2)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
                elseif ind(i) == 1
                    DIM = V.dim([2 3]);
                    C=[0 0 1 slice(1);1 0 0 0;0 1 0 0;0 0 0 1]; % Sagital Transformation Matrix
                    Xlim = [plotopts.lims.Ylim(1):plotopts.lims.Ylim(2)]; % Sagital Limits
                    Ylim = [plotopts.lims.Zlim(1):plotopts.lims.Zlim(2)];
                    nametotal = [ATtile '-Sagital_' num2str(slice(1))];
                    DTtitle = [ATtile ' : Sagital Slice Number ' num2str(slice(1)) '. Limits: Y from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
                end
                % Reading Background Image
                Is = spm_slice_vol(V,C,DIM,0);
                Is = Is(Xlim,Ylim);
                Is=rot90(Is);
                
                
                % Reading Overlay Image
                IOst = spm_slice_vol(VO(1),C,DIM,0);
                IOst = IOst(Xlim,Ylim);
                IOst = rot90(IOst);
                IOs(:,:,1) = IOst;
                
                IOst = spm_slice_vol(VO(2),C,DIM,0);
                IOst = IOst(Xlim,Ylim);
                IOst = rot90(IOst);
                IOs(:,:,2) = IOst;
                
                IOst = spm_slice_vol(VO(3),C,DIM,0);
                IOst = IOst(Xlim,Ylim);
                IOst = rot90(IOst);
                IOs(:,:,3) = IOst;
                
                % Plot Background
                imagesc(Is);
                colormap gray;axis equal;axis off;
                hold on;
                switch  plotopts.dtype.type
                    case 'scalar'
                        %% ============== Plotting Colored Maps ============ %%
                        h = imshow(IOs);
                        hold off;
                        set(h, 'AlphaData', plotopts.transp);
                    case 'vectorial'
                        %% ============ Plotting Vector Field ============== %%
                        indd = find((sum(IOs,3)~=0)&(~isnan(sum(IOs,3))));
                        [X,Y] = ind2sub([size(IOs,1) size(IOs,2)],indd);
                        It = squeeze(IOs(:,:,1));
                        drawev(:,1) = It(indd);
                        It = squeeze(IOs(:,:,2));
                        drawev(:,2) = It(indd);
                        It = squeeze(IOs(:,:,3));
                        drawev(:,3) = It(indd);
                        drawev(:,2) = -1*drawev(:,2);% 180 degrees rotation around X axis
                        if plotopts.dtype.modul == 0
                            norms = sqrt(sum(drawev.^2,2));
                            drawev = drawev./repmat(norms,[1 3]);
                        end
                        for k = 1:length(indd)
                            line([Y(k)-plotopts.linelength*drawev(k,1) Y(k)+plotopts.linelength*drawev(k,1)], ...
                                [X(k)-plotopts.linelength*drawev(k,2) X(k)+plotopts.linelength*drawev(k,2)], ...
                                [(X(k)*0+2)-plotopts.linelength*drawev(k,3) (X(k)*0+2)+plotopts.linelength*drawev(k,3)],...
                                'Color',abs(drawev(k,:)),'LineWidth',1);
                        end
                        clear drawev;
                end
            else
                %% ============ Plotting ODFs ========================== %%
                if ind(i) == 3
                    DIM = V.dim([1 2]);
                    C=[1 0 0 0;0 1 0 0;0 0 1 slice(3);0 0 0 1]; % Axial Transformation Matrix
                    Xlim = [plotopts.lims.Xlim(1):plotopts.lims.Xlim(2)]; % Axial Limits
                    Ylim = [plotopts.lims.Ylim(1):plotopts.lims.Ylim(2)];
                    nametotal = [ATtile '-Axial_' num2str(slice(3))];
                    DTtitle = [ATtile ' : Axial Slice Number ' num2str(slice(3)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Y from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
                elseif ind(i) == 2
                    DIM = V.dim([1 3]);
                    C=[1 0 0 0;0 0 1 slice(2);0 1 0 0;0 0 0 1]; % Coronal Transformation Matrix
                    Xlim = [plotopts.lims.Xlim(1):plotopts.lims.Xlim(2)]; % Coronal Limits
                    Ylim = [plotopts.lims.Zlim(1):plotopts.lims.Zlim(2)];
                    nametotal = [ATtile '-Coronal_' num2str(slice(2))];
                    DTtitle = [ATtile ' : Coronal Slice Number ' num2str(slice(2)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
                elseif ind(i) == 1
                    DIM = V.dim([2 3]);
                    C=[0 0 1 slice(1);1 0 0 0;0 1 0 0;0 0 0 1]; % Sagital Transformation Matrix
                    Xlim = [plotopts.lims.Ylim(1):plotopts.lims.Ylim(2)]; % Sagital Limits
                    Ylim = [plotopts.lims.Zlim(1):plotopts.lims.Zlim(2)];
                    nametotal = [ATtile '_Sag_' num2str(slice(1))];
                    DTtitle = [ATtile ' : Sagital Slice Number ' num2str(slice(1)) '. Limits: Y from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
                end
                
                % Reading Background Image
                Is = spm_slice_vol(V,C,DIM,0);
                Is = Is(Xlim,Ylim);
                Is=rot90(Is);
                imagesc(Is);
                colormap gray;axis equal;axis off;
                hold on;
                for k = 1:length(VO)
                    IOst = spm_slice_vol(VO(k),C,DIM,0);
                    indd = isnan(IOst);
                    IOst(indd) = 0;
                    IOst = IOst(Xlim,Ylim);
                    IOs(:,:,k) = rot90(IOst);
                end
                indd = find((sum(IOs,3)~=0)&(sum(IOs,3)~=NaN));
                [X,Y] = ind2sub([size(IOs,1) size(IOs,2)],indd);
                % ============ Creating Patch ============================%
                verts = [0 0 0]; F = [0 0 0]; Vt = [0 0 0];
                for k = 1:length(indd)
                    ODF = IOs(X(k),Y(k),:);ODF = ODF(:);
                    Vs = (repmat((1/2)*ODF./max(ODF),[1 3]).*dirs);
                    Vt = [Vt;Vs];
                    verts = [verts;Vs+ repmat( [Y(k) X(k) 2],[size(dirs,1) 1])];
                    F = [F;faces+max(F(:))];
                end
                Vt(1,:) = [];verts(1,:) = []; F(1,:) = [];
                Surf.SurfData.vertices = verts;
                Surf.SurfData.faces = F;
                a = Vt ./(repmat(sqrt(dot(Vt,Vt,2)),1,3) + eps);
                Surf.SurfData.FaceVertexCData =abs(a);
                patch('Vertices', Surf.SurfData.vertices, 'Faces', Surf.SurfData.faces,'FaceVertexCData', Surf.SurfData.FaceVertexCData,'FaceColor', 'interp', 'FaceLighting','phong',...
                    'LineStyle','-',...
                    'LineWidth',.05,...
                    'EdgeColor',[.3 .3 .3],...
                    'AmbientStrength',.4,...
                    'FaceLighting','phong',...
                    'SpecularColorReflectance',.2,...
                    'DiffuseStrength',.5,...
                    'BackFaceLighting', 'reverselit');
                camlight;
                axis equal
                %% ========== End of Plotting ODFs ======================%%
            end
            clear IOs;
            %% ============== End of Working With 4D Maps=============== %%
        end
        %% =============== Saving Figures ============================== %%
        title(DTtitle);
        Figurename = [OutDir filesep nametotal '.jpg'];
        print( gcf,'-djpeg','-r300',Figurename);
        a = getframe(hf);
        imwrite(a.cdata,Figurename,'JPEG')
        OutFigures  = strvcat(OutFigures,Figurename);
        close(hf);
    end
end
% ------------------ Removing Temporary Files ----------------------------%
% delete(deblank(OverlayImage(1,:)));
% delete(deblank(OverlayImage(1,:)));
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

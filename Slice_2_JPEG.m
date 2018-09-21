function OutFigures = Slice_2_JPEG(InputImage, slices, Outfile, plotopts);
%
% Syntax :
% OutFigures = Slice_2_JPEG(InputImage, slices, plotopts, Outfile);
%
% The following program shows how to view brain slices and to overlay maps
% on top of anatomical scans.
%
% Input Parameters:
%     InputImage       : Input Anatomical Image
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
%         Outfile        : Output Image prefix 
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
%InputImage = '/media/Data/PROCESSING_RESULTS/Plan_U/5-freesurfer_processing/s19enzo_h5c31dcp-20131015/tmp/T1_Masked.nii';
slices = [0 130 0];

% OverlayImage = strvcat('/home/yaleman/Data/7-connectome/ASPER_00001__101-20060510/dtreconst/dtk/ASPER_00001__101-20060510_odf.nii',...
% '/home/yaleman/Data/7-connectome/ASPER_00001__101-20060510/dtreconst/dtk/ASPER_00001__101-20060510_odfdirs.txt',...
% '/home/yaleman/Data/7-connectome/ASPER_00001__101-20060510/dtreconst/dtk/ASPER_00001__101-20060510_odffaces.txt');


%% ================== Cheking Input Parameters ========================= %%
close all;
V = spm_vol(InputImage);
[pth,name,ext] = fileparts(V(1).fname);
if nargin < 4
    plotopts.transp = 0.4;
    plotopts.colormap = 'jet';
    plotopts.itype = 'image';
    plotopts.atype = 'aparc+aseg';
    plotopts.dtype.type = 'vectorial';
    plotopts.dtype.modul = 0;
    plotopts.linelength = 1/2;
    %     plotopts.lims.Xlim = [20 120];
    %     plotopts.lims.Ylim = [20 120];
    %     plotopts.lims.Zlim = [10 50];
    OutDir = pth;
end
if nargin < 3
    Outfile = [pth filesep name];
end


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
TitleF = ['Image: ' upper(name) ];
ATtile = [ upper(name)];
OutFigures = '';
for t = 1:size(slices,1)
    slice = slices(t,:);
    ind = find(slice);
    for i = 1:length(ind)
        colordef black;
        hf = figure('numbertitle','off','name',TitleF,'Color',[0 0 0],'Position',[0 0 1200 900]);
        %% Non-Coloured Images
        if ind(i) == 3
            DIM = V(1).dim([1 2]);
            C=[1 0 0 0;0 1 0 0;0 0 1 slice(3);0 0 0 1]; % Axial Transformation Matrix
            Xlim = [plotopts.lims.Xlim(1):plotopts.lims.Xlim(2)]; % Axial Limits
            Ylim = [plotopts.lims.Ylim(1):plotopts.lims.Ylim(2)];
            nametotal = ['-Axial_' num2str(slice(3))];
            DTtitle = [ATtile ' : Axial Slice Number ' num2str(slice(3)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Y from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
        elseif ind(i) == 2
            DIM = V(1).dim([1 3]);
            C=[1 0 0 0;0 0 1 slice(2);0 1 0 0;0 0 0 1]; % Coronal Transformation Matrix
            Xlim = [plotopts.lims.Xlim(1):plotopts.lims.Xlim(2)]; % Coronal Limits
            Ylim = [plotopts.lims.Zlim(1):plotopts.lims.Zlim(2)];
            nametotal = ['-Coronal_' num2str(slice(2))];
            DTtitle = [ATtile ' : Coronal Slice Number ' num2str(slice(2)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
        elseif ind(i) == 1
            DIM = V(1).dim([2 3]);
            C=[0 0 1 slice(1);1 0 0 0;0 1 0 0;0 0 0 1]; % Sagital Transformation Matrix
            Xlim = [plotopts.lims.Ylim(1):plotopts.lims.Ylim(2)]; % Sagital Limits
            Ylim = [plotopts.lims.Zlim(1):plotopts.lims.Zlim(2)];
            nametotal = ['-Sagital_' num2str(slice(1))];
            DTtitle = [ATtile ' : Sagital Slice Number ' num2str(slice(1)) '. Limits: Y from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
        end
        % Reading Background Image
        Is = spm_slice_vol(V(1),C,DIM,0);
        Is = Is(Xlim,Ylim);% Selecting Limits
        Is=rot90(Is);
        
        % Plot Background
        imagesc(Is);colormap gray;axis equal;axis off;
        
        Title(DTtitle);
        Figurename = [Outfile nametotal '.jpg'];
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
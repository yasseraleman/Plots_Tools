function varargout = Prepare_ODFs_Slice_To_meshviewer(odfImage,slices,varargin);
%% ====================== Checking input parameters ===================== %
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return;
else
    % Reading ODF File or Variable
    temp = whos('odfImage');
    switch temp.class
        case 'char'
            if ~exist(odfImage,'file')
                error('The ODF file does not exist');
                return
            else
                [pthodf,nmodf,extodf] = fileparts(odfImage);
                switch deblank(extodf)
                    case '.gz'
                        boolzipODF = 1;
                        tempODF = unzip_nifti(odfImage);
                    case '.nii'
                        boolzipODF = 0;
                        tempODF = odfImage;
                    otherwise
                        error('Unrecognized image format.');
                        return;
                end
                Vodf = spm_vol(tempODF);
                
                param.Xlim = [1 Vodf(1).dim(1)];
                param.Ylim = [1 Vodf(1).dim(2)];
                param.Zlim = [1 Vodf(1).dim(3)];
                voxsize = sqrt(sum(Vodf(1).mat(1:3,1:3).^2));
            end
        case 'double'
            boolzipODF = 0;
    end
    
    % Images
    bgImage  = '';
    maskImage  = '';
    
    % Booleans
    boolzipMask = 0;
    boolzipMask = 0;
    
    
    % Parameters
    BGcolMap = 'gray'; % ColorMap
    
    
    % Recosntruction Scheme.
    filename = which('recon_matrix.m');
    [dirTemp, nameTemp] = fileparts(filename);
    [inputDir, diffToolsDir] = fileparts(dirTemp);
    param.reconsDir = [inputDir  filesep 'sampling_and_reconstruction_schemes' filesep 'On_the_sphere' filesep '724_shell.txt'];
    param.reconsFac = [inputDir  filesep 'sampling_and_reconstruction_schemes' filesep 'On_the_sphere' filesep '724_sphere_facets.txt'];
    
end


% deal with the input arguments
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return;
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'bgImage' % Background Image
                    bgImage=varargin{2};
                case 'maskImage' % Mask Image
                    maskImage=varargin{2};
                    
                    % Reading Mask Image or Variable
                    temp = whos('maskImage');
                    switch temp.class
                        case 'char'
                            if ~exist(maskImage,'file')
                                maskImage = '';
                            else
                                [pthmask,nmmask,extmask] = fileparts(maskImage);
                                switch deblank(extmask)
                                    case '.gz'
                                        boolzipMask = 1;
                                        tempMask = unzip_nifti(maskImage);
                                        Vmask = spm_vol(tempMask);
                                    case '.nii'
                                        boolzipMask = 0;
                                        tempMask = maskImage;
                                        Vmask = spm_vol(tempMask);
                                    otherwise
                                        warning('Unrecognized image format for Mask Image.');
                                        maskImage = '';
                                end
                            end
                        case 'double'
                            
                    end
                case 'xLims'
                    % Limits in X axis (Left-Right axis)
                    param.Xlim=varargin{2};
                case 'yLims'
                    % Limits in X axis (Anterior-Posterior axis)
                    param.Ylim=varargin{2};
                case 'zLims'
                    % Limits in Z axis (Inferior-Superior axis)
                    param.Zlim=varargin{2};
                case 'BGcolMap' % Background Image Colormap
                    BGcolMap=varargin{2};
                case 'reconsDir' % Reconstruction Scheme(Directions)
                    param.reconsDir=varargin{2};
                    if exist(param.reconsDir,'file')
                        dirs = load(deblank(param.reconsDir));
                        faces = convhulln(dirs);
                    else
                        warning('The reconstruction directions file does not exist. A default reconstruction scheme will be used.');
                    end
                case 'reconsFac' % Reconstruction Scheme(Faces)
                    param.reconsFac=varargin{2};
                    
                    if exist(param.reconsFac,'file')
                        faces = load(deblank(param.reconsFac));
                    else
                        warning('The reconstruction faces file does not exist. A default faces file will be used.');
                    end
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end
%% =================== End of checking input parameters ================= %

dirs = load(deblank(param.reconsDir));
faces = load(deblank(param.reconsFac));

% Creating Distance Neighborhood
gapSlicOdf = 0.5; 
compPeaks = 1;
Beta = 15;
norma = sqrt(sum(dirs.^2,2));
dirs = dirs./repmat(norma,[1 3]);
c = sqrt(2*(1-cos(Beta*pi/180)));
DistanceMatrix = dist(dirs');
ind = find(DistanceMatrix > c);
DistanceMatrix(ind) = 0;
Ndir = size(dirs,1);

allXs = [0;0];
allYs = [0;0];
allZs= [0;0];
%% ============================ Main Program ============================ %
cont = 0;
for t = 1:size(slices,1)
    slice = slices(t,:);
    ind = find(slice);
    for i = 1:length(ind)
        cont = cont + 1;
        
        %% ============ Plotting ODFs ========================== %%
        if ind(i) == 3
            DIM = Vodf(1).dim([1 2]);
            C=[1 0 0 0;0 1 0 0;0 0 1 slice(3);0 0 0 1]; % Axial Transformation Matrix
            Xlim = [param.Xlim(1):param.Xlim(2)]; % Axial Limits
            Ylim = [param.Ylim(1):param.Ylim(2)];
            %              nametotal = [ATtile '-Axial_' num2str(slice(3))];
            %             DTtitle = [ATtile ' : Axial Slice Number ' num2str(slice(3)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Y from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
        elseif ind(i) == 2
            DIM = Vodf(1).dim([1 3]);
            C=[1 0 0 0;0 0 1 slice(2);0 1 0 0;0 0 0 1]; % Coronal Transformation Matrix
            Xlim = [param.Xlim(1):param.Xlim(2)]; % Coronal Limits
            Ylim = [param.Zlim(1):param.Zlim(2)];
            %              nametotal = [ATtile '-Coronal_' num2str(slice(2))];
            %             DTtitle = [ATtile ' : Coronal Slice Number ' num2str(slice(2)) '. Limits: X from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
        elseif ind(i) == 1
            DIM = Vodf(1).dim([2 3]);
            C=[0 0 1 slice(1);1 0 0 0;0 1 0 0;0 0 0 1]; % Sagital Transformation Matrix
            Xlim = [param.Ylim(1):param.Ylim(2)]; % Sagital Limits
            Ylim = [param.Zlim(1):param.Zlim(2)];
            %              nametotal = [ATtile '_Sag_' num2str(slice(1))];
            %             DTtitle = [ATtile ' : Sagital Slice Number ' num2str(slice(1)) '. Limits: Y from ' num2str(min(Xlim)) ' to ' num2str(max(Xlim)) ' Z from ' num2str(min(Ylim)) ' to ' num2str(max(Ylim)) ];
        end
        
        % Reading ODF Image
        if exist('Vodf','var')
            for k = 1:length(Vodf)
                IOst = spm_slice_vol(Vodf(k),C,DIM,0);
                indd = isnan(IOst);
                IOst(indd) = 0;
                IOst = IOst(Xlim,Ylim);
                Iodf(:,:,k) = IOst;
                %Iodf(:,:,k) = rot90(IOst);
            end
        else
        end
        
        % Reading Mask Image
        if ~isempty(maskImage)
            
            if exist('Vmask','var')
                Imask = spm_slice_vol(Vmask(1),C,DIM,0);
                Imask = logical(Imask(Xlim,Ylim));
                %Imask = rot90(Imask);
            else
                temp = sum(Iodf,3);
                temp(isnan(temp)) = 0;
                Imask = logical(temp);
                
            end
        else
            temp = sum(Iodf,3);
            temp(isnan(temp)) = 0;
            Imask = logical(temp);
        end
        
        indd = find(Imask~=0);
        [X,Y] = ind2sub([size(Iodf,1) size(Iodf,2)],indd);
        
        
        %% ================= Creating ODF Patch ======================== %%
        ODF = zeros(size(dirs,1),length(indd));
        for k = 1:size(dirs,1)
            tempSlice = Iodf(:,:,k);
            ODF(k,:) = tempSlice(indd)';
        end
        ODF = ODF./repmat(max(ODF),[size(dirs,1) 1]);
        
        
        
        %% =================== Creating ODF Peaks ====================== %%
        if compPeaks
            PeaksMat = zeros(size(ODF));
            for z = 1:size(dirs,1)
                Neigh = unique(nonzeros(find(DistanceMatrix(z,:))));
                tempODF = repmat(ODF(z,:),[length(Neigh) 1]);
                DiffMat = tempODF - ODF(Neigh,:);
                
                index  = find(DiffMat <= 0); % Neighbors Vertices with smaller values
                DiffMat(index) = 0;
                DiffMat(isnan(DiffMat)) = 0;  % Removing Nan Values
                
                % Detecting Peaks
                DiffMat = sum(logical(DiffMat),1);
                index  = find(DiffMat ~=  length(Neigh));
                DiffMat(index) = 0;
                PeaksMat(z,:) = 0.5*logical(DiffMat).*ODF(z,:);
            end
            
            
            indpeaks = find(PeaksMat);
            Np = length(indpeaks);
            vertpeaks = zeros(length(indpeaks),6);

            if ind(i) == 3
                tempVar = repmat(X(:)',[size(dirs,1) 1]);tempVar = tempVar(:);
                vertpeaks(:,1) = tempVar(indpeaks);
                
                tempVar = repmat(Y(:)',[size(dirs,1) 1]);;tempVar = tempVar(:);
                vertpeaks(:,2) = tempVar(indpeaks);
                
                vertpeaks(:,3) = .5*ones(length(indpeaks),1) + slice(ind(i));
                
                tempVar = PeaksMat.*repmat(dirs(:,1),[1 length(indd)]) + repmat(X(:)',[size(dirs,1) 1]);tempVar = tempVar(:);
                vertpeaks(:,4) = tempVar(indpeaks);
                
                tempVar = PeaksMat.*repmat(dirs(:,2),[1 length(indd)]) + repmat(Y(:)',[size(dirs,1) 1]);tempVar = tempVar(:);
                vertpeaks(:,5) = tempVar(indpeaks);
                
                tempVar = PeaksMat.*repmat(dirs(:,3),[1 length(indd)]) + gapSlicOdf + slice(ind(i))*ones(size(dirs,1),length(X));tempVar = tempVar(:);
                vertpeaks(:,6) = tempVar(indpeaks);
                
            elseif ind(i) == 2
                tempVar = repmat(X(:)',[size(dirs,1) 1]);tempVar = tempVar(:);
                vertpeaks(:,1) = tempVar(indpeaks);
                
                vertpeaks(:,2) = .5*ones(length(indpeaks),1) + slice(ind(i));
                
                tempVar = repmat(Y(:)',[size(dirs,1) 1]);tempVar = tempVar(:);
                vertpeaks(:,3) = tempVar(indpeaks);
                
                tempVar = PeaksMat.*repmat(dirs(:,1),[1 length(indd)]) + repmat(X(:)',[size(dirs,1) 1]);tempVar = tempVar(:);
                vertpeaks(:,4) = tempVar(indpeaks);

                tempVar = PeaksMat.*repmat(dirs(:,2),[1 length(indd)]) + gapSlicOdf + slice(ind(i))*ones(size(dirs,1),length(X));tempVar = tempVar(:);
                vertpeaks(:,5) = tempVar(indpeaks);
                                
                tempVar = PeaksMat.*repmat(dirs(:,3),[1 length(indd)]) + repmat(Y(:)',[size(dirs,1) 1]);tempVar = tempVar(:);
                vertpeaks(:,6) = tempVar(indpeaks);
            elseif ind(i) == 1
                vertpeaks(:,1) = .5*ones(length(indpeaks),1)  + slice(ind(i));
                
                tempVar = repmat(X(:)',[size(dirs,1) 1]);tempVar = tempVar(:);
                vertpeaks(:,2) = tempVar(indpeaks);
                
                tempVar = repmat(Y(:)',[size(dirs,1) 1]);tempVar = tempVar(:);
                vertpeaks(:,3) = tempVar(indpeaks);
                
                tempVar = PeaksMat.*repmat(dirs(:,1),[1 length(indd)]) + gapSlicOdf + slice(ind(i))*ones(size(dirs,1),length(X));tempVar = tempVar(:);
                vertpeaks(:,4) = tempVar(indpeaks);
                
                tempVar = PeaksMat.*repmat(dirs(:,2),[1 length(indd)]) + repmat(X(:)',[size(dirs,1) 1]);tempVar = tempVar(:);
                vertpeaks(:,5) = tempVar(indpeaks);

                tempVar = PeaksMat.*repmat(dirs(:,3),[1 length(indd)]) + repmat(Y(:)',[size(dirs,1) 1]);tempVar = tempVar(:);
                vertpeaks(:,6) = tempVar(indpeaks);
            end
            
            allXs = [allXs vertpeaks(:,[1 4])'];
            allYs = [allYs vertpeaks(:,[2 5])'];
            allZs = [allZs vertpeaks(:,[3 6])'];
        end
        
        if ind(i) == 3
            Xdir = 0.4*ODF.*repmat(dirs(:,1),[1 length(indd)]) + repmat(X(:)',[size(dirs,1) 1]);
            Ydir = 0.4*ODF.*repmat(dirs(:,2),[1 length(indd)]) + repmat(Y(:)',[size(dirs,1) 1]);
            Zdir = 0.4*ODF.*repmat(dirs(:,3),[1 length(indd)]) + gapSlicOdf + slice(ind(i))*ones(size(dirs,1),length(X));
        elseif ind(i) == 2
            Xdir = 0.4*ODF.*repmat(dirs(:,1),[1 length(indd)]) + repmat(X(:)',[size(dirs,1) 1]);
            Ydir = 0.4*ODF.*repmat(dirs(:,2),[1 length(indd)]) + gapSlicOdf + slice(ind(i))*ones(size(dirs,1),length(X));
            Zdir = 0.4*ODF.*repmat(dirs(:,3),[1 length(indd)]) + repmat(Y(:)',[size(dirs,1) 1]);
        elseif ind(i) == 1
            Xdir = 0.4*ODF.*repmat(dirs(:,1),[1 length(indd)]) + gapSlicOdf + slice(ind(i))*ones(size(dirs,1),length(X));
            Ydir = 0.4*ODF.*repmat(dirs(:,2),[1 length(indd)]) + repmat(X(:)',[size(dirs,1) 1]);
            Zdir = 0.4*ODF.*repmat(dirs(:,3),[1 length(indd)]) + repmat(Y(:)',[size(dirs,1) 1]);
        end
        
        
% % % %         Xdir = 0.4*ODF.*repmat(dirs(:,1),[1 length(indd)]) + repmat(X(:)',[size(dirs,1) 1]);
% % % %         Ydir = 0.4*ODF.*repmat(dirs(:,2),[1 length(indd)]) + repmat(Y(:)',[size(dirs,1) 1]);
% % % %         Zdir = 0.4*ODF.*repmat(dirs(:,3),[1 length(indd)]) + gapSlicOdf + slice(ind(i))*ones(size(dirs,1),length(X));
        
        
        % Reordering ODF
        allIndexesODF = repmat([0:max(faces(:)):max(faces(:))*(length(indd)-1)],[size(faces,1) 1]); % Image indexes
        Xfac = allIndexesODF + repmat(faces(:,1),[1 length(indd)]);
        Yfac = allIndexesODF + repmat(faces(:,2),[1 length(indd)]);
        Zfac = allIndexesODF + repmat(faces(:,3),[1 length(indd)]);
        
        Surf.SurfData.vertices = [Xdir(:) Ydir(:) Zdir(:)]; clear Xdir Ydir Zdir;
        Surf.SurfData.faces    = [Xfac(:) Yfac(:) Zfac(:)]; clear Xfac Yfac Zfac;
        
        
        Xdir = ODF.*repmat(dirs(:,1),[1 length(indd)]);
        Ydir = ODF.*repmat(dirs(:,2),[1 length(indd)]);
        Zdir = ODF.*repmat(dirs(:,3),[1 length(indd)]);
        
        Surf.SurfData.FaceVertexCData  = abs([Xdir(:) Ydir(:) Zdir(:)]./(repmat(sqrt(dot([Xdir(:) Ydir(:) Zdir(:)],[Xdir(:) Ydir(:) Zdir(:)],2)),1,3) + eps));
        
        SurfODFs(cont) = Surf;
        %% ================= End of Creating Patch ===================== %%
    end
end

if length(SurfODFs) > 1
    for i = 1:length(SurfODFs);
        if i == 1
            Surfj.SurfData.vertices = SurfODFs(i).SurfData.vertices;
            Surfj.SurfData.faces = SurfODFs(i).SurfData.faces;
            Surfj.SurfData.FaceVertexCData = SurfODFs(i).SurfData.FaceVertexCData;
            
        else
            Surfj.SurfData.vertices = [Surfj.SurfData.vertices;SurfODFs(i).SurfData.vertices];
            Surfj.SurfData.faces =    [Surfj.SurfData.faces;SurfODFs(i).SurfData.faces+max(Surfj.SurfData.faces(:))];
            Surfj.SurfData.FaceVertexCData = [Surfj.SurfData.FaceVertexCData; SurfODFs(i).SurfData.FaceVertexCData];
        end
    end
    clear SurfODFs;
    SurfODFs = Surfj;
    clear Surfj;
end


%% ====================== Compressing Images =========================== %%
if  (boolzipODF == 1)&(exist('Vodf','var'))
    zip_nifti(Vodf(1).fname);
    delete(Vodf(1).fname);
    remove_niimat(Vodf(1).fname);
end
if  (boolzipMask == 1)&(exist('Vmask','var'))
    zip_nifti(Vmask(1).fname);
    delete(Vmask(1).fname);
    remove_niimat(Vmask(1).fname);
end
%% ====================== End of Compressing Images ==================== %%

%% ====================== Creating and Adding Background Image =========================== %%
if ~isempty(bgImage)
    Surfslices = Slice2Surf(bgImage,slices,'colMap',BGcolMap);
    
    for i = 1:length(Surfslices);
        SurfODFs.SurfData.vertices = [SurfODFs.SurfData.vertices;Surfslices(i).SurfData.vertices(:,2) Surfslices(i).SurfData.vertices(:,1) Surfslices(i).SurfData.vertices(:,3)];
        SurfODFs.SurfData.faces =    [SurfODFs.SurfData.faces;Surfslices(i).SurfData.faces+max(SurfODFs.SurfData.faces(:))];
        SurfODFs.SurfData.FaceVertexCData = [SurfODFs.SurfData.FaceVertexCData; Surfslices(i).SurfData.FaceVertexCData];
    end
end
%% =================== End of Creating and Adding Background Image ==================== %%
% Outputs
varargout{1} = SurfODFs;
if compPeaks
    allXs(:,1) = [];
    allYs(:,1) = [];
    allZs(:,1) = [];
    varargout{2} = [allXs(1,:)' allYs(1,:)' allZs(1,:)' allXs(2,:)' allYs(2,:)' allZs(2,:)'];
else
    varargout{2} = '';
end



return;

function varargout = Slice2Surf(bgImage,slices,varargin);
%
% Syntax :
% Surfout = Slice2Surf(bgImage,slices,colMap);
%
% This function creates a surface from a specified image slices.
%
% Input Parameters:
%   bgImage     : Image filename.
%   slices      : Slices (ie [0 0 128], [128 128 0])
%   colMap      : Colormap used to see the results.
%
% Output Parameters:
%   Surfout     : Surface of the slices.
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Red_Surf Plot_oversurf Atlas_Surf
% Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM
% March 30th 2017
% Version $1.0

%% ====================== Checking input parameters ===================== %
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return;
else
    % Reading ODF File or Variable
    temp = whos('bgImage');
    switch temp.class
        case 'char'
            if ~exist(bgImage,'file')
                error('The Image file does not exist');
                return
            else
                [pthbg,nmbg,extbg] = fileparts(bgImage);
                switch deblank(extbg)
                    case '.gz'
                        boolzipBG = 1;
                        tempBG = unzip_nifti(bgImage);
                    case '.nii'
                        boolzipBG = 0;
                        tempBG = bgImage;
                    otherwise
                        error('Unrecognized image format.');
                        return;
                end
                Vbg = spm_vol(tempBG);
                
                param.Xlim = [1 Vbg(1).dim(1)];
                param.Ylim = [1 Vbg(1).dim(2)];
                param.Zlim = [1 Vbg(1).dim(3)];
                voxsize = sqrt(sum(Vbg(1).mat(1:3,1:3).^2));
            end
        case 'double'
            boolzipBG = 0;
    end
    
    % Parameters
    colMap = 'jet'; % ColorMap
        
end


% deal with the input arguments
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return;
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'colMap' % Colormap
                    colMap=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end
%% =================== End of checking input parameters ================= %


%% ============================ Main Program ============================ %
cont = 0;
for t = 1:size(slices,1)
    slice = slices(t,:);
    ind = find(slice);
    for i = 1:length(ind)
        cont = cont + 1;
        
        switch temp.class
            case {'double','single'}
                
                if ind(i) == 3
                    Ibg = squeeze(bgImage(:,:,slice(3)));
                    sliceVal =  slice(3)-1;
                elseif ind(i) == 2
                    Ibg = squeeze(bgImage(:,slice(2),:));
                    sliceVal =  slice(2)-1;
                elseif ind(i) == 1
                    Ibg = squeeze(bgImage(slice(1),:,:));
                    sliceVal =  slice(1)-1;
                end
            otherwise
                %% ============ Detecting Slice ========================== %%
                if ind(i) == 3
                    DIM = Vbg(1).dim([1 2]);
                    C=[1 0 0 0;0 1 0 0;0 0 1 slice(3);0 0 0 1]; % Axial Transformation Matrix
                    Xlim = [param.Xlim(1):param.Xlim(2)]; % Axial Limits
                    Ylim = [param.Ylim(1):param.Ylim(2)];
                    sliceVal =  slice(3)-1;
                elseif ind(i) == 2
                    DIM = Vbg(1).dim([1 3]);
                    C=[1 0 0 0;0 0 1 slice(2);0 1 0 0;0 0 0 1]; % Coronal Transformation Matrix
                    Xlim = [param.Xlim(1):param.Xlim(2)]; % Coronal Limits
                    Ylim = [param.Zlim(1):param.Zlim(2)];
                    sliceVal =  slice(2)-1;
                elseif ind(i) == 1
                    DIM = Vbg(1).dim([2 3]);
                    C=[0 0 1 slice(1);1 0 0 0;0 1 0 0;0 0 0 1]; % Sagital Transformation Matrix
                    Xlim = [param.Ylim(1):param.Ylim(2)]; % Sagital Limits
                    Ylim = [param.Zlim(1):param.Zlim(2)];
                    sliceVal =  slice(1)-1;
                end
                
                %% ================== Creating Slice Surface ================== %%
                
                Ibg = spm_slice_vol(Vbg,C,DIM,0);
                Ibg = Ibg(Xlim,Ylim);
        end

        %Ibg=rot90(Ibg);
        dims = size(Ibg);
        [X,Y] = meshgrid([1:dims(1)],[1:dims(2)]);
        
        if ind(i) == 3
            Vertices = [X(:) Y(:) X(:)*0+sliceVal];
            Colors = Val2colors(Ibg(:),colMap);
        elseif ind(i) == 2
            Vertices = [X(:)*0+sliceVal X(:) Y(:)];
            Ibg = Ibg';
            Colors = Val2colors(Ibg(:),colMap);
        elseif ind(i) == 1
            Vertices = [X(:)  X(:)*0+sliceVal Y(:)];
            Ibg = Ibg';
            Colors = Val2colors(Ibg(:),colMap);
        end
        tri = delaunay(X,Y);
        Surfout(cont).SurfData.vertices = Vertices;
        Surfout(cont).SurfData.faces = flipdim(tri,2);
        Surfout(cont).SurfData.FaceVertexCData = Colors;

    end
    %% ================ End of Creating Slice Surface ================== %%
    
end

cont = 0;
for i = 1:length(Surfout);
    cont = cont + 1;
    if cont == 1
        Surfj.SurfData.vertices = Surfout(i).SurfData.vertices;
        Surfj.SurfData.faces = Surfout(i).SurfData.faces;
        Surfj.SurfData.FaceVertexCData = Surfout(i).SurfData.FaceVertexCData;

    else
        Surfj.SurfData.vertices = [Surfj.SurfData.vertices;Surfout(i).SurfData.vertices];
        Surfj.SurfData.faces =    [Surfj.SurfData.faces;Surfout(i).SurfData.faces+max(Surfj.SurfData.faces(:))];
        Surfj.SurfData.FaceVertexCData = [Surfj.SurfData.FaceVertexCData; Surfout(i).SurfData.FaceVertexCData];
    end
end
clear Surfout;
Surfout = Surfj;
clear Surfj;

%% ====================== Compressing Images =========================== %%
if  (boolzipBG == 1)&(exist('Vbg','var'))
    zip_nifti(Vbg(1).fname);
    delete(Vbg(1).fname);
    remove_niimat(Vbg(1).fname);
end
%% ====================== End of Compressing Images ==================== %%

%% ======================== End of Main Program ============================ %

% Outputs
varargout{1} = Surfout;
return;


function varargout = unzip_nifti(varargin);
%
% Syntax :
%      zipFilename = unzip_nifti(zipFilename, outDir);
%
% This script decompress nifty files.
%
% Input Parameters:
%        zipFilename               : Image filename.
%        outDir                    : Output directory.
%
%
% Output Parameters:
%        uzipFilename               : Unzipped Image.
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Cuban Neuroscience Center
% March 22th 2007
% Version $1.0

%% ============================= Checking Inputs ======================= %%

if nargin == 1
    zipFilename = varargin{1};
    [outDir, nm, ext] = fileparts(zipFilename);
    boolsave = 0;
elseif nargin == 2
    zipFilename = varargin{1};
    outDir      = varargin{2};
    boolsave = 0;
elseif nargin == 3
    zipFilename = varargin{1};
    outDir       = varargin{2};
    boolsave     = varargin{3};
else
    error('Please check input parameters');
    return;
end
if ~exist(zipFilename,'file')
    error('The specified zip file does not exist');
    return;
else
    [pth, nm, ext] = fileparts(zipFilename);
end
if isnumeric(boolsave)
    boolsave = logical(boolsave);
else
    warning('The boolsave variable must be numeric. The original file will be deleted: boolsave = 0');
    boolsave = 0;
end
%% ======================= End of Checking Inputs ====================== %%

%% ============================ Main Program =========================== %%
switch ext
    case '.gz'
        tempName = [outDir filesep nm];
        [pth, nm, ext] = fileparts(tempName);
        indnii = strfind(nm,'.nii');
        if ~isempty(indnii)
            if indnii == length(nm)-3;
                nm = nm(1:indnii-1);
            end
        end
        gunzip(zipFilename, outDir)
        uzipFilename =    [outDir filesep nm '.nii'];
        
        if ~boolsave
            delete(zipFilename);
        end
        
        %  Output
        varargout{1} = uzipFilename;
    case '.nii'
        disp('This file is not a zipped file.')
        uzipFilename = zipFilename;
        
        %  Output
        varargout{1} = uzipFilename;
end
%% ====================== End of Main Program ========================== %%
return

function Colors = Val2colors(txt,cl);
%
% Syntax :
% Colors = Val2colors(txt,cl);
%
% This function creates a colormap table according to the values of the txt
% variable. The colormap can be also specified.
%
% Input Parameters:
%      txt        : Vector of characteristic values
%       cl        : Colormap
%
%
% Output Parameter:
%     Colors      : Color Table
%
% See also: Surf_Color
%__________________________________________________
% Authors: Yasser Aleman Gomez
%LIM
% February 1th 2013
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ==1
    cl = 'jet';
end
[ut,i1t,it2] = unique(txt);
[col] = colormaps_colors(cl,size(ut,1));
Colors = col(it2,:);
return

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

function varargout = zip_nifti(varargin);
%
% Syntax :
%      zipFilename = zip_nifti(uzipFilename, outDir, boolsave);
%
% This script decompress nifty files.
%
% Input Parameters:
%        uzipFilename              : Image filename.
%        outDir                    : Output directory.
%        boolsave                  : Boolean variable to remove unzipped
%                                    image (1: Keep unzipped file,
%                                    0(default) otherwise).
%
% Output Parameters:
%        zipFilename               : Zipped Image.
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Cuban Neuroscience Center
% March 22th 2007
% Version $1.0

%% ============================= Checking Inputs ======================= %%

if nargin == 1
    uzipFilename = varargin{1};
    [outDir, nm, ext] = fileparts(uzipFilename);
    boolsave = 0;
elseif nargin == 2
    uzipFilename = varargin{1};
    outDir       = varargin{2};
    boolsave = 0;
elseif nargin == 3
    uzipFilename = varargin{1};
    outDir       = varargin{2};
    boolsave     = varargin{3};
else
    error('Please check input parameters');
    return;
end
if ~exist(uzipFilename,'file')
    error('The specified zip file does not exist');
    return;
else
    [pth, nm, ext] = fileparts(uzipFilename);
end
if isnumeric(boolsave)
    boolsave = logical(boolsave);
else
    warning('The boolsave variable must be numeric. The original file will be deleted: boolsave = 0');
    boolsave = 0;
end
%% ======================= End of Checking Inputs ====================== %%

%% ============================ Main Program =========================== %%
switch deblank(ext)
    case '.nii'
        gzip(uzipFilename, outDir);
        uzipFilename =    [outDir filesep nm '.nii'];
        
        %  Output
        varargout{1} = uzipFilename;
        if ~boolsave
            delete(uzipFilename);
        end
    case '.gz'
        disp('This file was previously zipped.')
        zipFilename = uzipFilename;
        
        %  Output
        varargout{1} = zipFilename;
end

%% ====================== End of Main Program ========================== %%
return

function remove_niimat(varargin);
%
% Syntax :
%      remove_niimat(niiFname);
%
% This script removes the mat file created  during the nifti saving
% process.
%
% Input Parameters:
%        niiFname                  : Image filename or Image volume (see spm_vol);
%
% Output Parameters:
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Cuban Neuroscience Center
% February 21th 2006
% Version $1.0

niiFname = varargin{1};
[pth, nm, ext] = fileparts(niiFname);

if exist([pth filesep nm filesep '.mat']);
    delete([pth filesep nm filesep '.mat']);
end
return;

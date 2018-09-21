function varargout = remove_outliers(varargin);
%
% Syntax :
%     dataOut = remove_outliers(dataIn);
%
% This function removes outliers from a 3D data matrix. 
% 3D matrix are asumed as Nstruc x Nmethods x Nsubjects
% 2D matrix are asumed as Nsubjects x Nstruc or Nmethods
%
% Input Parameters:
%        dataIn                         : Input Data
%        opts (optional input)          : Options:
%                                         opts.noutliers - Number of outliers to be removed
%                                         opts.confinterv - Confidence interval (0-1)
%                                         
%
% Output Parameters:
%         dataOut                       : Output Data
%
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% November 13th 2014
% Version $1.0

%% ============================= Checking Inputs ======================= %%
if nargin == 0
    error('Two Inputs are mandatory');
    return
elseif nargin == 1
    dataIn = varargin{1};
    opts.noutliers = 3;
    opts.confinterv = 0.95;
elseif nargin == 2
    dataIn = varargin{1};
    opts = varargin{2};
elseif nargin >2
    error('To Many Input Parameters');
    return;
end

if nargout > 1
    error('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%

%% =========================== Main Program ============================ %%
dataSize = size(dataIn);
dimData = length(dataSize);
dataOut = zeros(dataSize);
if dimData == 3
    for i = 1:dataSize(1)
        for j = 1:dataSize(2)
            X = squeeze(dataIn(i,j,:));
            Y = squeeze(dataIn(i,j,:));
            
            indinf = find(X == Inf );
            indnan = find(isnan(X));
            indrealx = find(X ~=Inf);
            indrealx(find(ismember(indrealx,indnan)==1))=[];
            X(unique([indinf indnan])) =  mean(X(indrealx));
            
            indinf = find(Y == Inf );
            indnan = find(isnan(Y));
            indrealy = find(Y ~=Inf);
            indrealy(find(ismember(indrealy,indnan)==1))=[];
            Y(unique([indinf indnan])) =  mean(Y(indrealx));
            
            
            
            exy = confellipse2([X Y],opts.confinterv);
            in = inpolygon(X,Y,exy(:,1),exy(:,2));
            outliers_idx = find(in == 0);
            
            %[X,Y,rSquares,outliers_idx] = regoutliers(X,Y,noutliers,plotOp);
            %X(outliers_idx) =Inf;
            
            X(outliers_idx)  = mean(X(indrealx));
            if find(X > 40)
                a = 1;
            end
            dataOut(i,j,:) = X;
        end
    end
elseif dimData == 2
    if dataSize(1) == 1
        dataIn = dataIn(:);
    end
    for j = 1:dataSize(2)
        X = squeeze(dataIn(:,j));
        Y = squeeze(dataIn(:,j));
        
        indinf = find(isfinite(X) ==0);
        indrealx = find(isfinite(X) ==1);
        X(indinf) =  mean(X(indrealx));
        indinf = find(isfinite(Y) ==0);
        indreal = find(isfinite(Y) ==1);
        Y(indinf) =  mean(Y(indreal));
        
        
        exy = confellipse2([X Y],.95);
        in = inpolygon(X,Y,exy(:,1),exy(:,2));
        outliers_idx = find(in == 0);
        
        %[X,Y,rSquares,outliers_idx] = regoutliers(X,Y,noutliers,plotOp);
        X(outliers_idx) =Inf;
        X(outliers_idx)  = mean(X(indrealx));
        dataOut(:,j) = X;
    end
else
    dataOut = dataIn;
end
%% ========================= End of Main Program ======================= %%
% Outputs
varargout{1} = dataOut;
return;
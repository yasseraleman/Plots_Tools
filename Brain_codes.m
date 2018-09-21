function [GMcodes,Names, Colors] = Brain_codes(Image, atlastype);
%
% Syntax :
% [GMcodes,Names, Colors] = Brain_codes(atlastype);
%
% This function extract all gray matter codes from a specified atlas.
%
% Input Parameters:
%   atlastype     : Atlas type.
%
%
% Output Parameters:
%   GMcodes       : Gray matter codes.
%   Names         : Structure Names
%   Colors        : Structure Colors
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
% March 20th 2012
% Version $1.0
GMcodes = '';
Names = '';
cont = 0;
Colors = [0 0 0];
V0 = spm_vol(Image);
IOs = spm_read_vols(V0);
GMcodes = unique(IOs(:));
GMcodes(GMcodes==0) = [];
switch atlastype
    case 'aparc+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        %GMcodes = [7 8 10:13 17:18  26 46 47 49:54 58 1001:1003 1005:1035 2001:2003 2005:2035];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'a2009s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        %GMcodes = [7 8 10:13 17:18  26 46 47 49:54 58 11101:11175 12101:12175];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'a2005s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        %GMcodes = [7 8 10:13 17:18  26 46 47 49:54 58 1102:1181 2102:2181];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'yeoh7networks'
        
    case 'yeoh17networks'
        
    case 'parckmeans'
        
     case 'ibaspm116'
        txtfile = which('atlas116.cod');
        col = 255*[1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
        Ncolor = size(col,1);
        
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        re = floor(length(Codes)/Ncolor); col = repmat(col,[re+1 1]);
        %         Codes = str2num(char(a));
        Names = char(b);
        Colors = col(1:length(Codes),:);
        GMcodes = Codes;
    case 'ibaspm71'
        txtfile = which('atlas71.cod');
        col = 255*[1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
        Ncolor = size(col,1);
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        %         Codes = str2num(char(a));
        Names = char(b);
        re = floor(length(Codes)/Ncolor); col = repmat(col,[re+1 1]);
        Colors = col(1:length(Codes),:);
        GMcodes = Codes;
    otherwise
        
end

return
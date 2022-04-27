% /* --------------------------------------------------------------------------------------
%  * File:    MigDescriptor2movit.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

% Create file to visualize descriptors in BioEmergences Movit

% Descriptors Migration
% DescriptorsMigration=[cellid cellnum speed speedAve speedAveDirection_v neigh negihchange cellDiam Convergence]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Dataset parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadParametersKinMovit;%import parameters to run the translation to movit

%Load the raw descriptors
load(rawVelDescriptors);
fStep=min( DescriptorsV(:,3));
lStep=max( DescriptorsV(:,3));

size(DescriptorsV)

%indexes to convert to colormaps and all that
sp=7;
spmean=11
spr=15
neigh=24;%each eigenvalue are 5 values-> norm, sign, vector
neighGain=25;
cellDiam=26;%each eigenvalue are 5 values-> norm, sign, vector
convergence=27;
sl=19
slmeanX=20
slmeanY=21
slmeanZ=22
slmeanM=23
indexToNormalize=[sp spmean spr sl slmeanX slmeanY slmeanZ slmeanM neigh neighGain cellDiam convergence]

isSymmetric=[0 0 0 0 0 0 0 0 0 0 0 0]%0-no symmetric %1-symmetric by max %2- symmetric by min
startInZero=[1 1 1 1 0 0 0 1 1 1 1 1]
mxper=[90 90 90 90 90 90 90 90 90 90 90 90]
mnper=[0 0 0 0 0 0 0 0 0 0 0 0]
ms=255;
mis=0;


%Calculate the global for all timesteps max and min values, we use a percentil for it.
existent=find( DescriptorsV(:,sp)>-1);
Whole= DescriptorsV(existent,:);
size(Whole)

%Find the global max and min for the whole sequence!!
%We need to skip using -1
pmatrix=Whole(:,indexToNormalize);
%values have been discarded already
[maxp minp]= getMaxMinDesc(pmatrix, mxper, mnper, -100);

%Create the header passing the global max and min.
if vheader==1
    fid=fopen(fileHeader,'w');
    %writeHeaderMigration_v1(fid, temp_ave, spat_ave, maxp, minp, bins);
    %fclose(fid);
elseif vheader==2
    %Different way having local max and min for all the steps
    fid2=fopen(fileHmax,'w');
    fid3=fopen(fileHmin,'w');
    fid=fopen(fileHeader,'w');
    %writeHeaderMigration_v12(fid, temp_ave, spat_ave, maxp, minp, bins)
elseif vheader==3
    %writeHeader();%to design and implement
end
%fclose(fid);

allmax=zeros((tfinal-tini+1), size(indexToNormalize,2));
allmin=zeros((tfinal-tini+1), size(indexToNormalize,2));
cm=1;
for i=tini:tfinal
    stamp=i+1
    fileoutstep=[tlapse 't' num2str(i) '.kdr'];
    fid1=fopen(fileoutstep, 'w');
    it=find( DescriptorsV(:,3)==(i)); %we stored as zero-based like in movit... should be save in zero based too, but next version
    Ts= DescriptorsV(it,:);
    sT=size(Ts)
    existent=find( Ts(:,sp)>-1 );%~((T(:,pp)>-1) && (T(:,pn)>-1)));
    %size(existent)
    Timestep = Ts(existent,:);
    
    tpmatrix=Timestep(:,indexToNormalize);
    [maxpt minpt]= getMaxMinDesc(tpmatrix, mxper, mnper, -100);
    allmax(cm,:)=maxpt;
    allmin(cm,:)=minpt;
    cm=cm+1;
    
    if vheader==2
        %We save the max/min of the timestep in a file next to the global
        %header for vheader==2.
        fprintf(fid2, '%s',['-1;-1;-1;' num2str(maxpt(1)) ';-1;-1;-1;' num2str(maxpt(2)) ';-1;-1;-1;' num2str(maxpt(3)) ';-1;-1;-1;' ...
            num2str(maxpt(4)) ';' num2str(maxpt(5)) ';' num2str(maxpt(6)) ';' num2str(maxpt(7)) ';' num2str(maxpt(8)) ...
            ';' num2str(maxpt(9)) ';' num2str(maxpt(10)) ';' num2str(maxpt(11)) ';' num2str(maxpt(12))]);
        fprintf(fid2, '\n');
        fprintf(fid3, '%s',['-1;-1;-1;' num2str(minpt(1)) ';-1;-1;-1;' num2str(minpt(2)) ';-1;-1;-1;' num2str(minpt(3)) ';-1;-1;-1;' ...
            num2str(minpt(4)) ';' num2str(minpt(5)) ';' num2str(minpt(6)) ';' num2str(minpt(7)) ';' num2str(minpt(8)) ...
            ';' num2str(minpt(9)) ';' num2str(minpt(10)) ';' num2str(minpt(11)) ';' num2str(minpt(12))]);
        fprintf(fid3, '\n');
    end
    
    %Max and minimum of the timestep regarding the global max and min
    %We do it by index!!
    if indexNorm
        imaxt=getColorIndexDesc(maxpt, maxp, minp, bins, -100);
        imint=getColorIndexDesc(minpt, maxp, minp, bins, -100);
    else
        imaxt=maxpt;
        imint=minpt;
    end
    
    fprintf(fid1, '%s', ['cellid;cellnum;timestep;maxspeed;maxspeedave;maxspeedavex;maxspeedavey;maxspeedavez;maxsl;maxmeanslx;maxmeansly;maxmeanslz;maxmeanslm;maxeneigh;maxneighgain;maxeneighX;maxneighgainX;']);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s',['-1;-1;-1;' num2str(imaxt(1)) ';-1;-1;-1;' num2str(imaxt(2)) ';-1;-1;-1;' num2str(imaxt(3)) ';-1;-1;-1;' ...
        num2str(imaxt(4)) ';' num2str(imaxt(5)) ';' num2str(imaxt(6)) ';' num2str(imaxt(7)) ';' num2str(imaxt(8)) ...
        ';' num2str(imaxt(9)) ';' num2str(imaxt(10)) ';' num2str(imaxt(11)) ';' num2str(imaxt(12))]);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s', ['cellid;cellnum;timestep;minspeed;minspeedave;minspeedavex;minspeedavey;minspeedavez;minsl;minmeanslx;minmeansly;minmeanslz;minmeanslm;mineneigh;minneighgain;mineneighX;minneighgainX']);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s',['-1;-1;-1;' num2str(imint(1)) ';-1;-1;-1;' num2str(imint(2)) ';-1;-1;-1;' num2str(imint(3)) ';-1;-1;-1;' ...
        num2str(imint(4)) ';' num2str(imint(5)) ';' num2str(imint(6)) ';' num2str(imint(7)) ';' num2str(imint(8)) ...
        ';' num2str(imint(9)) ';' num2str(imint(10)) ';' num2str(imint(11)) ';' num2str(imint(12))]);
    fprintf(fid1, '\n');
    
    %Normalization with the global max and min physical vales
    %Then we convert to index within the colormap!
    for ic=1:size(Timestep,1)
        if globalNorm
            iCol=getColorIndexDesc(Timestep(ic,indexToNormalize), maxp, minp, bins, -100);
        else
            iCol=getColorIndexDesc(Timestep(ic,indexToNormalize), maxpt, minpt, bins, -100);
        end
        cellid=Timestep(ic,1);
        cellnum=Timestep(ic,2);
        % Timestep(ic,e+1)
        %  Timestep(ic,c+2)
        %   Timestep(ic,rota+3)
        
        fprintf(fid1, '%s', [num2str(cellid) ';' num2str(cellnum) ';' num2str(i) ';' num2str(iCol(1)) ';' ...
            num2str(Timestep(ic,sp+1)) ';' num2str(Timestep(ic,sp+2)) ';' num2str(Timestep(ic,sp+3)) ';' ...
            num2str(iCol(2)) ';' num2str(Timestep(ic,spmean+1)) ';' num2str(Timestep(ic,spmean+2)) ';' num2str(Timestep(ic,spmean+3)) ';' ...
            num2str(iCol(3)) ';' num2str(Timestep(ic,spr+1)) ';' num2str(Timestep(ic,spr+2)) ';' num2str(Timestep(ic,spr+3)) ';' ...
            num2str(iCol(4)) ';' num2str(iCol(5)) ';' num2str(iCol(6)) ';' num2str(iCol(7)) ';' num2str(iCol(8)) ';' ...
            num2str(iCol(9)) ';' num2str(iCol(10)) ';' num2str(iCol(11)) ';' num2str(iCol(12)) ';']);
        fprintf(fid1, '\n');
        
    end   
end

maxp=max(allmax);
minp=min(allmin);
writeHeaderMigration_v12(fid, temp_ave, spat_ave, maxp, minp, bins);
%fclose(fid);
fclose(fid1);

if vheader==2
    fclose(fid2);
    fclose(fid3);
end

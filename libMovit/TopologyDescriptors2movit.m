% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

% Create file to visualize descriptors in BioEmergences Movit

% Descriptors Topology to Movit
% speed 
% speedAve 
% Topology 
% P
% Q 
% RotationDiscriminat 
% Expansion 
% Compression
% Rotation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DescriptorsTopology=[cellid cellnum t x y z speed speedAve topology P Q RotationDiscriminat Expansion exp_nv Compression comp_nv Rotation rota_nv]

%Load the raw descriptors
load(rawTopoDescriptors);
%load(rawVelDescriptors);
fStep=min( DescriptorsT(:,3));
lStep=max( DescriptorsT(:,3));
%indexes to convert to colormaps and all that
speed=7; speedAve=8; q=14; topo=12;
p=13;rotDisc=15;e=16;c=20; rota=24;s=topo;

%Normalize
indexToNormalize=[speed speedAve p q rotDisc e c rota];
isSymmetric=[0 0 1 0 0 0 0 0]%0-no symmetric %1-symmetric by max %2- symmetric by min
startInZero=[1 1 0 0 0 1 1 1]
mxper=[90 90 85 85 85 90 90 90]
mnper=[0 0 15 15 10 0 0 0 ]
ms=255;
mis=0;

%Create the header passing the global max and min.
if vheader==1
    fid=fopen(fileHeader,'w');
    % writeHeaderTopo_v12(fid, temp_ave, spat_ave, maxp, minp, bins);
    %fclose(fid);
elseif vheader==2
    %Different way having local max and min for all the steps
    fid2=fopen(fileHmax,'w');
    fid3=fopen(fileHmin,'w');
    fid=fopen(fileHeader,'w');
    % writeHeaderTopo_v12(fid, temp_ave, spat_ave, maxp, minp, bins)
elseif vheader==3
    % writeHeader();%to design and implement
end
%fclose(fid);

%Save local max and min
allmax=zeros((tfinal-tini+1), size(indexToNormalize,2));
allmin=zeros((tfinal-tini+1), size(indexToNormalize,2));
cm=1;

for i=tini:tfinal
    ['Exporting Topology Descriptors timestep ' num2str(i)]
    fileoutstep=[tlapse 't' num2str(i) '.kdr'];
    fid1=fopen(fileoutstep, 'w');
    it=find( DescriptorsT(:,3)==(i));
    Ts= DescriptorsT(it,:);
    sT=size(Ts)
    mst=ms;
    mist=mis;
    Timestep=Ts;
    clear Ts;
    
    tpmatrix=Timestep(:,indexToNormalize);
    [maxpt minpt]= getMaxMinDesc(tpmatrix, mxper, mnper, -1);%dont consider the -1
    %symmetric
    %maxpt(4)
    %minpt(4)
    maxpt(isSymmetric==1)= max(abs(maxpt(isSymmetric==1)),abs(minpt(isSymmetric==1)));
    minpt(isSymmetric==1)= -max(abs(maxpt(isSymmetric==1)),abs(minpt(isSymmetric==1)));
    maxpt(isSymmetric==2)= min(abs(maxpt(isSymmetric==2)),abs(minpt(isSymmetric==2)));
    minpt(isSymmetric==2)= -min(abs(maxpt(isSymmetric==2)),abs(minpt(isSymmetric==2)));
    minpt(startInZero==1)= 0;
    %save
    allmax(cm,:)=maxpt;
    allmin(cm,:)=minpt;
   
    cm=cm+1;

    if vheader==2
        %We save the max/min of the timestep in a file next to the global
        %header for vheader==2.
        fprintf(fid2, '%s',['-1;-1;' num2str(i) ';' num2str(maxpt(1))  ';' num2str(maxpt(2)) ';-1;-1;-1;' num2str(mst) ';' ...
            num2str(maxpt(3)) ';' num2str(maxpt(4)) ';' ...
            num2str(maxpt(5)) ';' num2str(maxpt(6)) ';-1;-1;-1;'  num2str(maxpt(7)) ';-1;-1;-1;' num2str(maxpt(8)) ';-1;-1;-1']);
        fprintf(fid2, '\n');
        fprintf(fid3, '%s',['-1;-1;' num2str(i) ';' num2str(minpt(1)) ';' num2str(minpt(2)) ';-1;-1;-1;' num2str(mist) ';' ...
            num2str(minpt(3)) ';' num2str(minpt(4)) ';' ...
            num2str(minpt(5)) ';' num2str(minpt(6)) ';-1;-1;-1;'  num2str(minpt(7)) ';-1;-1;-1;' num2str(minpt(8)) ';-1;-1;-1' ]);
        fprintf(fid3, '\n');
    end
    
    %Max and minimum of the timestep regarding the global max and min
    if indexNorm
        imst=getColorIndexDesc(mst, ms, mis, bins);
        imist=getColorIndexDesc(mist, ms, mis, bins);
        imaxt=getColorIndexDesc(maxpt, maxp, minp, bins);
        imint=getColorIndexDesc(minpt, maxp, minp, bins);
        if StartOne
            imaxt=imaxt+1;
            imint=imint+1;
        end
    else
        imst=mst;
        imist=mist;
        imaxt=maxpt;
        imint=minpt;
    end
    
    fprintf(fid1, '%s', ['cellid;cellnum;timestep;maxspeed;maxspeedave;maxspeedx;maxspeedy;maxspeedz;maxtopo;maxp;maxq;maxD;maxexp;maxexpx;maxexpy;maxexpz;maxcomp;maxcompx;maxcompy;maxcompz;maxrot;maxrotx;maxroty;maxrotz']);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s',['-1;-1;-1;' num2str(imaxt(1)) ';' num2str(imaxt(2)) ';-1;-1;-1;' num2str(imst) ';' num2str(imaxt(3)) ...
        ';' num2str(imaxt(4)) ';' num2str(imaxt(5)) ';' num2str(imaxt(6)) ';-1;-1;-1;'  ...
        num2str(imaxt(7)) ';-1;-1;-1;' num2str(imaxt(8)) ';-1;-1;-1;' ] );
    fprintf(fid1, '\n');
    fprintf(fid1, '%s', ['cellid;cellnum;timestep;minspeed;minspeedave;minspeedx;minspeedy;minspeedz;mintopo;minq;minp;minD;minexp;minexpx;minexpy;minexpz;mincomp;mincompx;mincompy;mincompz;minrot;minrotx;minroty;minrotz']);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s',['-1;-1;-1;' num2str(imint(1)) ';' num2str(imint(2)) ';-1;-1;-1;'  num2str(imist) ';' num2str(imint(3)) ...
        ';' num2str(imint(4)) ';' num2str(imint(5)) ';' num2str(imint(6)) ';-1;-1;-1;'  ...
        num2str(imint(7)) ';-1;-1;-1;' num2str(imint(8)) ';-1;-1;-1;'] );
    fprintf(fid1, '\n');
    
    %Normalization with the global max and min physical vales
    %Then we convert to index within the colormap!
    for ic=1:size(Timestep,1)
        if globalNorm
            is=getColorIndexDesc(Timestep(ic,s), ms, mis, bins);
            iCol=getColorIndexDesc(Timestep(ic,indexToNormalize), maxp, minp, bins);
        else
            is=getColorIndexDesc(Timestep(ic,s), mst, mist, bins);
            iCol=getColorIndexDesc(Timestep(ic,indexToNormalize), maxpt, minpt, bins);
        end
        if StartOne
           iCol=iCol+1;
        end
        cellid=Timestep(ic,1);
        cellnum=Timestep(ic,2);
       
        fprintf(fid1, '%s', [num2str(cellid) ';' num2str(cellnum) ';' num2str(i) ';' num2str(iCol(1)) ';' ...
            num2str(iCol(2)) ';' num2str(Timestep(ic,speedAve+1)) ';' num2str(Timestep(ic,speedAve+2)) ';' num2str(Timestep(ic,speedAve+3)) ...
            ';' num2str(is) ';' num2str(iCol(3)) ';' num2str(iCol(4)) ';' num2str(iCol(5)) ...
            ';' num2str(iCol(6)) ';' num2str(Timestep(ic,e+1)) ';' num2str(Timestep(ic,e+2)) ';' num2str(Timestep(ic,e+3)) ';'  ...
            num2str(iCol(7)) ';' num2str(Timestep(ic,c+1)) ';' num2str(Timestep(ic,c+2)) ';' num2str(Timestep(ic,c+3)) ';' ...
            num2str(iCol(8)) ';' num2str(Timestep(ic,rota+1)) ';' num2str(Timestep(ic,rota+2)) ';' num2str(Timestep(ic,rota+3))]);
        fprintf(fid1, '\n');
    end
end

maxp=max(allmax);
minp=min(allmin);
   
%maxp(4)
%minp(4)
maxp(isSymmetric==1)= max(abs(maxp(isSymmetric==1)),abs(minp(isSymmetric==1)));
minp(isSymmetric==1)= -max(abs(maxp(isSymmetric==1)),abs(minp(isSymmetric==1)));
maxp(isSymmetric==2)= min(abs(maxp(isSymmetric==2)),abs(minp(isSymmetric==2)));
minp(isSymmetric==2)= -min(abs(maxp(isSymmetric==2)),abs(minp(isSymmetric==2)));
minp(startInZero==1)= 0;
 
writeHeaderTopo_v12(fid, temp_ave, spat_ave, maxp, minp, bins);
fclose(fid1);
if vheader==2
    fclose(fid2);
    fclose(fid3);
end


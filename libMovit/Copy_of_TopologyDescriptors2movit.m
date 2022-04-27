% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

% Create file to visualize descriptors in BioEmergences Movit

% Descriptors Topology to movit
% DescriptorsTopology=[cellid cellnum t x y z speed speedAve neigh s p q e exp_nv c comp_nv rota rota_nv]


%Load the raw descriptors
load(rawTopoDescriptors);
%load(rawVelDescriptors);
fStep=min( DescriptorsT(:,3));
lStep=max( DescriptorsT(:,3));
%indexes to convert to colormaps and all that
sp=7; spAve=8; neigh=12; s=13;
p=14;q=15;e=16;c=20;rota=24;
indexToNormalize=[sp spAve neigh p q e c rota]
%Calculate the global for all timesteps max and min values, we use a percentil for it.
ms=255;
mis=0;
%existent=find( DescriptorsT(:,s)>-1);%~((DescriptorsT(:,pp)>-1) && (DescriptorsT(:,pn)>-1)) );
%Whole= DescriptorsT(existent,:);
%we need the speed even if the F is not evaluated.
Whole=DescriptorsT;
size(Whole)
DescriptorsT;

%Find the global max and min
% pmatrix=Whole(:,indexToNormalize);
% size(pmatrix)
% [maxp minp]= getMaxMinDesc(pmatrix, mxper, mnper, -1);

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

allmax=zeros((tfinal-tini+1), size(indexToNormalize,2));
allmin=zeros((tfinal-tini+1), size(indexToNormalize,2));
cm=1;
for i=tini:tfinal
    stamp=i+1
    fileoutstep=[tlapse 't' num2str(i) '.kdr'];
    fid1=fopen(fileoutstep, 'w');
    it=find( DescriptorsT(:,3)==(i));
    %we stored as zero-based like in movit... should be save in zero based too, but next version
    Ts= DescriptorsT(it,:);
    sT=size(Ts)
    %We get the timestep max and min for each descriptor
    mst=ms;
    mist=mis;
    %existent=find( Ts(:,s)>-1 );%~((T(:,pp)>-1) && (T(:,pn)>-1)));
    %Timestep = Ts(existent,:);
    Timestep=Ts;
    clear Ts;
    
    tpmatrix=Timestep(:,indexToNormalize);
    [maxpt minpt]= getMaxMinDesc(tpmatrix, mxper, mnper, -1);
    allmax(cm,:)=maxpt;
    allmin(cm,:)=minpt;
    cm=cm+1;
    %       sp=7; spAve=8; neigh=12; s=13;
    %       p=14;q=15;e=16;c=20;rota=24;
    if vheader==2
        %We save the max/min of the timestep in a file next to the global
        %header for vheader==2.
        fprintf(fid2, '%s',['-1;-1;' num2str(i) ';' num2str(maxpt(1))  ';' num2str(maxpt(2)) ';-1;-1;-1;' num2str(maxpt(3)) ';' num2str(mst) ';' num2str(maxpt(4)) ';' num2str(maxpt(5)) ';' ...
            num2str(maxpt(6)) ';-1;-1;-1;'  num2str(maxpt(7)) ';-1;-1;-1;' num2str(maxpt(8)) ';-1;-1;-1']);
        fprintf(fid2, '\n');
        fprintf(fid3, '%s',['-1;-1;' num2str(i) ';' num2str(minpt(1)) ';' num2str(minpt(2)) ';-1;-1;-1;' num2str(minpt(3)) ';' num2str(mist) ';' num2str(minpt(4)) ';' ...
            num2str(minpt(5)) ';' num2str(minpt(6)) ';-1;-1;-1;'  num2str(minpt(7)) ';-1;-1;-1;' num2str(minpt(8)) ';-1;-1;-1' ]);
        fprintf(fid3, '\n');
    end
    
    %Max and minimum of the timestep regarding the global max and min
    %We do it by index!!
    if indexNorm
        imst=getColorIndexDesc(mst, ms, mis, bins);
        imist=getColorIndexDesc(mist, ms, mis, bins);
        imaxt=getColorIndexDesc(maxpt, maxp, minp, bins);
        imint=getColorIndexDesc(minpt, maxp, minp, bins);
    else
        imst=mst;
        imist=mist;
        imaxt=maxpt;
        imint=minpt;
    end
    
    fprintf(fid1, '%s', ['cellid;cellnum;timestep;maxspeed;maxspeedave;maxspeedx;maxspeedy;maxspeedz;maxq;maxsign;maxp;maxq;maxexp;maxexpx;maxexpy;maxexpz;maxcomp;maxcompx;maxcompy;maxcompz;maxrot;maxrotx;maxroty;maxrotz']);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s',['-1;-1;-1;' num2str(imaxt(1)) ';' num2str(imaxt(2)) ';-1;-1;-1;' num2str(imaxt(3)) ';' num2str(imst) ...
        ';' num2str(imaxt(4)) ';' num2str(imaxt(5)) ';' num2str(imaxt(6)) ';-1;-1;-1;'  ...
        num2str(imaxt(7)) ';-1;-1;-1;' num2str(imaxt(8)) ';-1;-1;-1;' ] );
    fprintf(fid1, '\n');
    fprintf(fid1, '%s', ['cellid;cellnum;timestep;minspeed;minspeedave;minspeedx;minspeedy;minspeedz;minq;minsign;minp;maxq;minexp;minexpx;minexpy;minexpz;mincomp;mincompx;mincompy;mincompz;minrot;minrotx;minroty;minrotz']);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s',['-1;-1;-1;' num2str(imint(1)) ';' num2str(imint(2)) ';-1;-1;-1;' num2str(imint(3)) ';' num2str(imist) ...
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
        cellid=Timestep(ic,1);
        cellnum=Timestep(ic,2);
        % Timestep(ic,e+1)
        %  Timestep(ic,c+2)
        %   Timestep(ic,rota+3)
        fprintf(fid1, '%s', [num2str(cellid) ';' num2str(cellnum) ';' num2str(i) ';' num2str(iCol(1)) ';' num2str(iCol(2)) ';' num2str(Timestep(ic,spAve+1)) ';' num2str(Timestep(ic,spAve+2)) ';' num2str(Timestep(ic,spAve+3)) ...
            ';' num2str(iCol(3)) ';' num2str(is) ';' num2str(iCol(4)) ';' num2str(iCol(5)) ...
            ';' num2str(iCol(6)) ';' num2str(Timestep(ic,e+1)) ';' num2str(Timestep(ic,e+2)) ';' num2str(Timestep(ic,e+3)) ';'  ...
            num2str(iCol(7)) ';' num2str(Timestep(ic,c+1)) ';' num2str(Timestep(ic,c+2)) ';' num2str(Timestep(ic,c+3)) ';' ...
            num2str(iCol(8)) ';' num2str(Timestep(ic,rota+1)) ';' num2str(Timestep(ic,rota+2)) ';' num2str(Timestep(ic,rota+3))]);
        fprintf(fid1, '\n');
        
    end
end
maxp=max(allmax);
minp=min(allmin);
writeHeaderTopo_v12(fid, temp_ave, spat_ave, maxp, minp, bins);
fclose(fid1);
if vheader==2
    fclose(fid2);
    fclose(fid3);
end


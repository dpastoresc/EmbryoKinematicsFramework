% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

% Create file to visualize descriptors in BioEmergences Movit

% Export Lagrangian Descriptors to movit
% Lagrange=[cellid cellnum t x y z FTLE isoFTLE FTLEback isoFTLEback Rotat e1 e2 e3 MaxShear]

%Load the raw descriptors
load(lDescriptors);

if left
    %fStep=min( DescriptorsLLeft(:,3))
    %lStep=max( DescriptorsLLeft(:,3))
    %DescriptorsLLeft(:,3)=50;
    fStep=min( DescriptorsLLeft(:,3));
    lStep=max( DescriptorsLLeft(:,3));
else
    fStep=min( DescriptorsL(:,3));
    lStep=max( DescriptorsL(:,3));
end
%indexes to convert to colormaps and all that
LCSex=7;    LCSdt=8;   backLCSex=9;    backLCSdt=10;
s=11;  volRatio=12;  rot=13; distor=29;
e1=17;  e2=21; e3=25;
ms=255;
mis=0;
indexToNormalize=[LCSex LCSdt backLCSex backLCSdt volRatio rot e1 e2 e3 distor];

if left
    Whole= DescriptorsLLeft;
    clear DescriptorsLLeft
else
    Whole= DescriptorsL;%we already removed all the nondesired stuff
    clear DescriptorsL
end
size(Whole)

%Find the global max and min for the whole sequence!! Polarity is
%preserved. No taking absolute magnitudes.
%pmatrix=Whole(:,indexToNormalize);
%ineg=pmatrix<0;    %pmatrix=abs(pmatrix);
%size(pmatrix)
%the bad rows are just not stored, so we dont need -1 anywhere
%[maxp minp]= getMaxMinDesc(pmatrix, mxper, mnper,-99999999);

%Create the header passing the global max and min.
if vheader==1
    fid=fopen(fileHeader,'w');
    writeHeaderLagrange_v12(fid, temp_ave, spat_ave, maxp, minp, bins, T);
    %fclose(fid);
elseif vheader==2
    %Different way having local max and min for all the steps
    fid2=fopen(fileHmax,'w');
    fid3=fopen(fileHmin,'w');
    fid=fopen(fileHeader,'w');
elseif vheader==3
    writeHeader();%to design and implement
end

tfinal2=tfinal-T+1
allmax=zeros((tfinal2-tini+1), size(indexToNormalize,2));
allmin=zeros((tfinal2-tini+1), size(indexToNormalize,2));
cm=1;
tini
tfinal
%pause

for i=tini:tfinal%(tfinal+T-1)%tini
    stamp=i+1;
    'Exporting lag to movit step'
    i
    it=find( Whole(:,3)==(i)); %we stored as zero-based like in movit... should be save in zero based too, but next version
    ceellsStep=size(it)
    %pause
    if size(it,1)>0
        
        mst=ms;
        mist=mis;
        fileoutstep=[tlapse 't' num2str(i) '.kdr'];
        fid1=fopen(fileoutstep, 'w');
        Ts= Whole(it,:);
        sT=size(Ts);
        %existent=find( T(:,S)>-1 );%~((T(:,pp)>-1) && (T(:,pn)>-1)));
        %Timestep = T(existent,:);
        Timestep=Ts;
        clear Ts;
        %We make the absolute to calculate values of the colormap
        %tpmatrix=abs(Timestep(:,indexToNormalize));
        tpmatrix=(Timestep(:,indexToNormalize));
        [maxpt minpt]= getMaxMinDesc(tpmatrix, mxper, mnper,-99999999);
        allmax(cm,:)=maxpt;
        allmin(cm,:)=minpt;
        cm=cm+1;
        
        size(maxpt)
        
        if vheader==2
            %We save the max/min of the timestep in a file next to the global
            %header for vheader==2.
            fprintf(fid2, '%s',['-1;-1;-1;' num2str(maxpt(1)) ';' num2str(maxpt(2)) ';' num2str(maxpt(3)) ';' num2str(maxpt(4)) ';'  num2str(mst) ';' ...
                num2str(maxpt(5)) ';' num2str(maxpt(6)) ';-1;-1;-1;'  num2str(maxpt(7)) ';-1;-1;-1;' ...
                num2str(maxpt(8)) ';-1;-1;-1;' num2str(maxpt(9)) ';-1;-1;-1;' num2str(maxpt(10))]);
            fprintf(fid2, '\n');
            fprintf(fid3, '%s',['-1;-1;-1;' num2str(minpt(1)) ';' num2str(minpt(2)) ';' num2str(minpt(3)) ';' num2str(minpt(4)) ';' num2str(mist) ';' ...
                num2str(minpt(5)) ';' num2str(minpt(6)) ';-1;-1;-1;'  num2str(minpt(7)) ';-1;-1;-1;' ...
                num2str(minpt(8)) ';-1;-1;-1;' num2str(minpt(9)) ';-1;-1;-1;' num2str(minpt(10))]);
            fprintf(fid3, '\n');
        end
        
        %Max and minimum of the timestep regarding the global max and min
        %We do it by index!!
        if indexNorm
            imst=getColorIndexDesc(mst, ms, mis, bins);
            imist=getColorIndexDesc(mist, ms, mis, bins);
            imaxt=getColorIndexDesc(minpt, maxp, minp, bins,-99999);
            imint=getColorIndexDesc(minpt, maxp, minp, bins,-99999);
        else
            imst=mst;
            imist=mist;
            imaxt=maxpt;
            imint=minpt;
        end
        
        fprintf(fid1, '%s', ['cid;cnum;tstep;maxLCS;maxisoLCS;maxVelo;maxs;maxexpansionRate;maxinv2;maxrotation_m;maxrotation_x;maxrotation_y;maxrotation_z;' ...
            'maxe1_m;maxe1_x;maxe1_y;maxe1_z;maxe2_m;maxe2_x;maxe2_y;maxe2_z;maxe3_m;maxe3_x;maxe3_y;maxe3_z;maxShear']);
        fprintf(fid1, '\n');
        fprintf(fid1, '%s',['-1;-1;-1;' num2str(imaxt(1)) ';' num2str(imaxt(2)) ';' num2str(imaxt(3)) ';' num2str(imaxt(4)) ';' num2str(imst) ...
            ';' num2str(imaxt(5)) ';' num2str(imaxt(6)) ';-1;-1;-1;' num2str(imaxt(7)) ';-1;-1;-1;' ...
            num2str(imaxt(8)) ';-1;-1;-1;' num2str(imaxt(9)) ';-1;-1;-1;' num2str(imaxt(10))]);
        fprintf(fid1, '\n');
        fprintf(fid1, '%s', ['cid;cnum;tstep;minLCS;minisoLCS;minVelo;mins;minexpansionRate;mininv2;minrotation_m;minrotation_x;minrotation_y;minrotation_z;' ...
            'mine1_m;mine1_x;mine1_y;mine1_z;mine2_m;mine2_x;mine2_y;mine2_z;mine3_m;mine3_x;mine3_y;mine3_z;minShear']);
        fprintf(fid1, '\n');
        fprintf(fid1, '%s',['-1;-1;-1;' num2str(imint(1)) ';' num2str(imint(2)) ';' num2str(imint(3)) ';' num2str(imint(4)) ';' num2str(imist) ...
            ';' num2str(imint(5)) ';' num2str(imint(6)) ';-1;-1;-1;' num2str(imint(7)) ';-1;-1;-1;' ...
            num2str(imint(8)) ';-1;-1;-1;' num2str(imint(9)) ';-1;-1;-1;' num2str(imint(10))]);
        fprintf(fid1, '\n');
        
        %Normalization with the global max and min physical vales
        %Then we convert to index within the colormap!
        for ic=1:size(Timestep,1)
            
            values=Timestep(ic,indexToNormalize);
            if globalNorm
                is=getColorIndexDesc(Timestep(ic,s), ms, mis, bins, -9999);
                iCol=getColorIndexDesc(values, maxp, minp, bins, -9999);
            else
                is=getColorIndexDesc(Timestep(ic,s), mst, mist, bins, -9999);
                iCol=getColorIndexDesc(values, maxpt, minpt, bins, -9999);
            end
            cellid=Timestep(ic,1);
            cellnum=Timestep(ic,2);
            
            %num2str(Timestep(ic,(e1+1)))
            %num2str(Timestep(ic,e1+1))
            
            
            fprintf(fid1, '%s', [num2str(cellid) ';' num2str(cellnum) ';' num2str(i) ';' num2str(iCol(1)) ';' num2str(iCol(2)) ...
                ';' num2str(iCol(3)) ';' num2str(iCol(4)) ';' num2str(is) ';' ...
                num2str(iCol(5)) ';' num2str(iCol(6)) ';' num2str(Timestep(ic,(rot+1))) ';' num2str(Timestep(ic,rot+2)) ';' num2str(Timestep(ic,rot+3)) ';' ...
                num2str(iCol(7)) ';' num2str(Timestep(ic,(e1+1))) ';' num2str(Timestep(ic,e1+2)) ';' num2str(Timestep(ic,e1+3)) ';' ...
                num2str(iCol(8)) ';' num2str(Timestep(ic,(e2+1))) ';' num2str(Timestep(ic,e2+2)) ';' num2str(Timestep(ic,e2+3)) ';' ...
                num2str(iCol(9)) ';'  num2str(Timestep(ic,(e3+1))) ';' num2str(Timestep(ic,e3+2)) ';' num2str(Timestep(ic,e3+3)) ';' num2str(iCol(10))]);
            fprintf(fid1, '\n');
            
        end
        fclose(fid1);
    end
end
maxp=max(allmax);
minp=min(allmin);
if vheader==2
    writeHeaderLagrange_v12(fid, temp_ave, spat_ave, maxp, minp, bins, T)
    fclose(fid2);
    fclose(fid3);
end


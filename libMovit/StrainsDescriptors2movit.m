% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

% Create file to visualize descriptors in BioEmergences Movit

% Descriptors Strains
% Qs 
% e1 e2 e3 
% Qd 
% Ratio=Qs-Qd 
% MaxShearAngle 
% de1 de2 de3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DescriptorsStrains=[cellid cellnum t x y z Qs e1 e2 e3 Qd Ratio MaxShearAngle de1 de2 de3]

%Load the raw descriptors
load(rawStrainDescriptors);
fStep=min( DescriptorsS(:,3));
lStep=max( DescriptorsS(:,3));
%indexes to convert to colormaps and all that
Qs=7; Qd=20; MS=22; R=21;
e1=8;e2=12;e3=16;d1=23;d2=27;d3=31;

%Normalize
indexToNormalize=[Qs e1 e2 e3 Qd R MS d1 d2 d3]
isSymmetric=[0 1 1 1 0 0 0 1 1 1]
startInZero=[1 0 0 0 1 0 1 0 0 0]
mxper=[85 90 90 90 85 90 90 90 90 90]
mnper=[0 10 10 10 0 10 0 10 10 10]

%Create the header passing the global max and min.
if vheader==1
    fid=fopen(fileHeader,'w');
    %writeHeaderStrains_v12(fid, temp_ave, spat_ave, maxp, minp, bins);
    %fclose(fid);
elseif vheader==2
    %Different way having local max and min for all the steps
    fid2=fopen(fileHmax,'w');
    fid3=fopen(fileHmin,'w');
    fid=fopen(fileHeader,'w');
    %writeHeaderStrains_v12(fid, temp_ave, spat_ave, maxp, minp, bins)
elseif vheader==3
    %writeHeader();%to design and implement
end
%fclose(fid);

%Save local max and min
allmax=zeros((tfinal-tini+1), size(indexToNormalize,2));
allmin=zeros((tfinal-tini+1), size(indexToNormalize,2));
cm=1;

for i=tini:tfinal
    ['Exporting Strain Rate Descriptors timestep ' num2str(i)]
    fileoutstep=[tlapse 't' num2str(i) '.kdr'];
    fid1=fopen(fileoutstep, 'w');
    it=find( DescriptorsS(:,3)==(i));
    Ts= DescriptorsS(it,:);
    sT=size(Ts)
    %Remove undesired rows first
    existent=find( Ts(:,Qs)>-1 );
    Timestep = Ts(existent,:);
    clear Ts;

    tpmatrix=Timestep(:,indexToNormalize);
    [maxpt minpt]= getMaxMinDesc(tpmatrix, mxper, mnper, -9999999999);
    %symmetric
    maxpt(isSymmetric==1)= max(abs(maxpt(isSymmetric==1)),abs(minpt(isSymmetric==1)));
    minpt(isSymmetric==1)= -max(abs(maxpt(isSymmetric==1)),abs(minpt(isSymmetric==1)));
    minpt(startInZero==1)= 0;
    %save
    allmax(cm,:)=maxpt;
    allmin(cm,:)=minpt;
    cm=cm+1;
    
    if vheader==2
        %We save the max/min of the timestep in a file next to the global
        %header for vheader==2.
        fprintf(fid2, '%s',['-1;-1;' num2str(i) ';' num2str(maxpt(1)) ';' ...
            num2str(maxpt(2)) ';-1;-1;-1;' num2str(maxpt(3)) ';-1;-1;-1;' ...
            num2str(maxpt(4)) ';-1;-1;-1;' num2str(maxpt(5)) ';'  num2str(maxpt(6)) ';' ...
            num2str(maxpt(7)) ';' num2str(maxpt(8)) ';-1;-1;-1;' num2str(maxpt(9)) ...
            ';-1;-1;-1;' num2str(maxpt(10)) ';-1;-1;-1;' ]);
        fprintf(fid2, '\n');
        fprintf(fid3, '%s',['-1;-1;' num2str(i) ';' num2str(minpt(1)) ';' ...
            num2str(minpt(2)) ';-1;-1;-1;' num2str(minpt(3)) ';-1;-1;-1;' ...
            num2str(minpt(4)) ';-1;-1;-1;' num2str(minpt(5)) ';'  num2str(minpt(6)) ';' ...
            num2str(minpt(7)) ';' num2str(minpt(8)) ';-1;-1;-1;' num2str(minpt(9)) ...
            ';-1;-1;-1;' num2str(minpt(10)) ';-1;-1;-1;' ]);
        fprintf(fid3, '\n');
    end
    
    %Max and minimum of the timestep regarding the global max and min
    %We do it by index!!
    if indexNorm
        imaxt=getColorIndexDesc(maxpt, maxp, minp, bins,-2);
        imint=getColorIndexDesc(minpt, maxp, minp, bins,-2);
        if StartOne
            imaxt=imaxt+1;
            imint=imint+1;
        end
    else
        imaxt=maxpt;
        imint=minpt;
    end
    
    fprintf(fid1, '%s', ['cellid;cellnum;timestep;maxs;maxe1;maxe1x;maxe1y;maxe1z;maxe2;maxe2x;maxe2y;maxe2z;maxe3;maxe3x;maxe3y;maxe3z;maxD;maxMS;maxR;maxe1;maxe1x;maxe1y;maxe1z;maxe2;maxe2x;maxe2y;maxe2z;maxe3;maxe3x;maxe3y;maxe3z']);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s',['-1;-1;-1;' num2str(imaxt(1)) ';' num2str(imaxt(2)) ';-1;-1;-1;' num2str(imaxt(3)) ';-1;-1;-1;' ...
        num2str(imaxt(4)) ';-1;-1;-1;' num2str(imaxt(5)) ';' num2str(imaxt(6)) ';' num2str(imaxt(7)) ';' ...
        num2str(imaxt(8)) ';-1;-1;-1;' num2str(imaxt(9)) ';-1;-1;-1;' num2str(imaxt(10)) ';-1;-1;-1']);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s', ['cellid;cellnum;timestep;mins;mine1;mine1x;mine1y;mine1z;mine2;mine2x;mine2y;mine2z;mine3;mine3x;mine3y;mine3z;minD;minMS;minR;minde1;minde1x;min1y;minde1z;minde2;minde2x;minde2y;minde2z;minde3;minde3x;minde3y;minde3z']);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s',['-1;-1;-1;' num2str(imint(1)) ';' num2str(imint(2)) ';-1;-1;-1;' num2str(imint(3)) ';-1;-1;-1;' ...
        num2str(imint(4)) ';-1;-1;-1;' num2str(imint(5)) ';' num2str(imint(6)) ';' num2str(imint(7)) ';' ...
        num2str(imint(8)) ';-1;-1;-1;' num2str(imint(9)) ';-1;-1;-1;' num2str(imint(10)) ';-1;-1;-1']);
    fprintf(fid1, '\n');

    %Normalization with the global max and min physical vales
    %Then we convert to index within the colormap!
    for ic=1:size(Timestep,1)
        
        values=Timestep(ic,indexToNormalize);
        if globalNorm
            %-2 should not affect as except Value
            iCol=getColorIndexDesc(values, maxp, minp, bins, -2);
        else
            iCol=getColorIndexDesc(values, maxpt, minpt, bins, -2);
        end
        if StartOne
            iCol=iCol+1;
        end
        cellid=Timestep(ic,1);
        cellnum=Timestep(ic,2);
        
        fprintf(fid1, '%s', [num2str(cellid) ';' num2str(cellnum) ';' num2str(i) ';' num2str(iCol(1)) ';' ...
            num2str(iCol(2)) ';' num2str(Timestep(ic,e1+1)) ';' num2str(Timestep(ic,(e1+2))) ';' num2str(Timestep(ic,(e1+3))) ';' ...
            num2str(iCol(3)) ';' num2str(Timestep(ic,(e2+1))) ';' num2str(Timestep(ic,e2+2)) ';' num2str(Timestep(ic,(e2+3))) ';' ...
            num2str(iCol(4)) ';' num2str(Timestep(ic,(e3+1))) ';' num2str(Timestep(ic,e3+2)) ';' num2str(Timestep(ic,e3+3)) ';' ...
            num2str(iCol(5)) ';' num2str(iCol(6)) ';' num2str(iCol(7)) ';' ...
            num2str(iCol(8)) ';' num2str(Timestep(ic,d1+1)) ';' num2str(Timestep(ic,d1+2)) ';' num2str(Timestep(ic,d1+3)) ';' ...
            num2str(iCol(9)) ';' num2str(Timestep(ic,d2+1)) ';' num2str(Timestep(ic,d2+2)) ';' num2str(Timestep(ic,d2+3)) ';' ...
            num2str(iCol(10)) ';' num2str(Timestep(ic,d3+1)) ';' num2str(Timestep(ic,d3+2)) ';' num2str(Timestep(ic,d3+3)) ]);
        fprintf(fid1, '\n');
        
    end
end

maxp=max(allmax);
minp=min(allmin);
maxp(isSymmetric==1)= max(abs(maxp(isSymmetric==1)),abs(minp(isSymmetric==1)));
minp(isSymmetric==1)= -max(abs(maxp(isSymmetric==1)),abs(minp(isSymmetric==1)));
minp(startInZero==1)= 0;
writeHeaderStrains_v12(fid, temp_ave, spat_ave, maxp, minp, bins);
fclose(fid1);
if vheader==2
    fclose(fid2);
    fclose(fid3);
end


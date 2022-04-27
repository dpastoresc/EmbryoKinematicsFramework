% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

% Create file to visualize descriptors in BioEmergences Movit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Dataset parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadParametersKinMovit;%import parameters to run the translation to movit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
descM=[desc 'M' seltag tagDesc '-' param];
mkdir(descM);
fileHeader=[descM '/header.kdh']
tlapse=[descM '/timelapse/'];
mkdir(tlapse);
if vheader==2
    fileHmax=[descM '/maxtemp.kdh']
    fileHmin=[descM '/mintemp.kdh']
end

%Load the raw descriptors
%load(rawMechanomeDescriptors);
fStep=min( DescriptorsM(:,3));
lStep=max( DescriptorsM(:,3));

size(DescriptorsM)
indexToNormalize=[7:size(DescriptorsM,2)]

%Calculate the global for all timesteps max and min values, we use a percentil for it.
existent=find( DescriptorsM(:,7)>-1);
Whole= DescriptorsM(existent,:);
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
    it=find( DescriptorsM(:,3)==(i)); %we stored as zero-based like in movit... should be save in zero based too, but next version
    Ts= DescriptorsM(it,:);
    sT=size(Ts)
    existent=find( Ts(:,7)>-1 );%~((T(:,pp)>-1) && (T(:,pn)>-1)));
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
        fprintf(fid2, '%s',['-1;-1;-1;' num2str(maxpt(1)) ';' num2str(maxpt(2)) ';' num2str(maxpt(3)) ...
            ';' num2str(maxpt(4)) ';' num2str(maxpt(5)) ';' num2str(maxpt(6))]);
        fprintf(fid2, '\n');
        fprintf(fid3, '%s',['-1;-1;-1;' num2str(minpt(1)) ';' num2str(minpt(2)) ';' num2str(minpt(3)) ...
            ';' num2str(minpt(4)) ';' num2str(minpt(5)) ';' num2str(minpt(6))]);
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
    
    fprintf(fid1, '%s', ['cellid;cellnum;timestep;maxmode1;maxmode2;maxmode3;maxmode4;maxmode5;maxmode6;']);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s',['-1;-1;-1;' num2str(imaxt(1)) ';' num2str(imaxt(2)) ';' num2str(imaxt(3)) ...
        ';' num2str(imaxt(4)) ';' num2str(imaxt(5)) ';' num2str(imaxt(6))]);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s', ['cellid;cellnum;timestep;minmode1;minmode2;minmode3;minmode4;minmode5;minmode6;']);
    fprintf(fid1, '\n');
    fprintf(fid1, '%s',['-1;-1;-1;' num2str(imint(1)) ';' num2str(imint(2)) ';' num2str(imint(3)) ...
        ';' num2str(imint(4)) ';' num2str(imint(5)) ';' num2str(imint(6))]);
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
        
        fprintf(fid1, '%s', [num2str(cellid) ';' num2str(cellnum) ';' num2str(i) ';' num2str(iCol(1)) ';' num2str(iCol(2)) ';' ...
             num2str(iCol(3)) ';' num2str(iCol(4)) ';' num2str(iCol(5)) ';' num2str(iCol(6))]);
        fprintf(fid1, '\n');
        
    end   
end

maxp=max(allmax);
minp=min(allmin);
writeHeaderMechanome_v12(fid, temp_ave, spat_ave, maxp, minp, bins);
%fclose(fid);
fclose(fid1);

if vheader==2
    fclose(fid2);
    fclose(fid3);
end

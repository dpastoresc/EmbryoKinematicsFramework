% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%%%%%%%%%%%%%%%%% SPATIO-TEMPORAL CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%

addpath('../colormaps')
addpath('colormaps')
color_cluster={'phase','phase','phase','phase','J','smart','phase','phase','phase'}
    
%%%% Load coherent trajectories
%track_file=[desc 'tracking-lag-' num2str(tinilag) '-' num2str(tfinal) '.csv']
track_file=[desc 'tracking-' num2str(tinilag) '-' num2str(tfinal) '.csv']
Track_lag=dlmread(track_file);
lag_steps=[max(t_start,tinilag):min(t_limit,tfinal)];%It starts in the points we want
mat_ids=[min(Track_lag(:,1)):max(Track_lag(:,1))];%The ids
mat_ids=[unique(Track_lag(:,1))];
num_ids=length(mat_ids)

%Careful, here we consider consistent ids! Otherwise we should select by
%time
t_ref=tinilag
if only_selection
    'Restrict Data to the Selections'
    %do the clustering only for the selections u want
    %if a temporal reference is given, the ids are selected in there
    if t_ref>-1
        selec_t=selec(selec(:,3)==t_ref,:);
        Track_lag_t=Track_lag(Track_lag(:,3)==t_ref,:);
        Desc_t=Desc(Desc(:,3)==t_ref,:);
    else
        Track_lag_t=Track_lag;
    end
    
    %see the spatial ids for the selectionss
    selec_i= arrayfun(@(x) find(selec_t(:,2)==x), clusterSelec, 'UniformOutput', false);
    size(selec_i)
    selec_ii=[]
    for i=1:length(selec_i)
        selec_ii=vertcat(selec_ii,selec_i{i});
    end
    
    %look in the tracking to get the corrsponing material_ids
    spatial_ids=selec_t(selec_ii,1);
    length(spatial_ids);
    spatial_ids=unique(spatial_ids)
    length(spatial_ids);
    %take the rows from Desc of the material ids within the selections
    selec2= arrayfun(@(x) find(Track_lag_t(:,2)==x), spatial_ids, 'UniformOutput', false);
    size(selec2);
    selec2= cell2mat(selec2);
    sTrack_lag=Track_lag_t(selec2,:);
else
    sTrack_lag=Track_lag;
end

material_ids=unique(sTrack_lag(:,1));
num_ids=length(material_ids)

%One dimension for the cluster
%mat_ids * steps * {time/spatial_id/desc}
ClusterData=zeros(num_ids,(lag_steps(end)-lag_steps(1)+1),3);
ClusterData(:)=NaN;%we want only complete trajectories

%Crop in time
Track_lag=Track_lag((Track_lag(:,3)>=t_start) & (Track_lag(:,3)<=t_limit),:);

%Go through all the material data to create the Cluster Datas
toremove=[];
ce=0;
for ii=1:num_ids
    
    'Material id profile'
    m_id=material_ids(ii)
    %take the trajectory->spatio-temporal data
    traj=Track_lag(Track_lag(:,1)==m_id,2:3);
    
    %Now we have the spatial ids for the material_id
    if DescriptorsSet>3
        %as this was built with the tracking, we just find for the m_id
        selecidsc= find(Desc(:,2)==m_id );%the cell num was substituted by the track_id  
        'Lag track' 
        size(selecidsc)
        selecids= selecidsc;        
        if length(selecids)~=size(ClusterData,2)
            ce=ce+1
            toremove(ce)=ii
            'ERROR'
        else
            data=Desc(selecids,[3 1 descriptor]);
            ClusterData(ii,:,:)=data;
        end
       % pause
    else
        selecidsc= arrayfun(@(x) find(Desc(:,1)==x), traj(:,1), 'UniformOutput', false);    
        size(selecidsc)
        selecids= cell2mat(selecidsc);
        size(selecids)
        if size(selecids,1)~=size(ClusterData,2)
            'ERROR'
            pause
        else
            data=Desc(selecids,[3 1 descriptor]);
            ClusterData(ii,:,:)=data;
        end
    end
end

ClusterData(toremove,:,:)=[];%remove strange things
if t_max_cluster>-1
    itcluster=find(ClusterData(1,:,1)==t_max_cluster)
else
    itcluster=100;
end
ClusterData=ClusterData(:,1:itcluster,:);

ce
num_ids
size(ClusterData)
min_t=min(ClusterData(1,:,1));
max_t=max(ClusterData(1,:,1));
time_axis=[min_t:10:max_t];
%Remove uncomplete trajectories (the time has been already cropped)
toremove=[];
c=1;
numOfOriginalIds=size(ClusterData,1)
for iii=1:numOfOriginalIds
    %lacks an id??
    if any(isnan(ClusterData(iii,:,2)))
        %remove iii from clus
        toremove(c)=iii;
        c=c+1;
    end
end

%Just take the 2D data for the clusering
ClusterData(toremove,:,:)=[];
clus=squeeze(ClusterData(:,:,3));%just the data but we keep ids

%Remove invalid outliars in a material perspective that will affect the clustering
%Filtering to smooth if necessary
%ave=[]
%cave=1
toremove=[];
c=1;
numOfMaterialProfiles=size(clus)
%general statistical cropping
mx=prctile(clus(:),98)
mn=prctile(clus(:),5)
clus(clus>mx)=mx;
clus(clus<mn)=mn;
smoothN=20
%     ave(cave)=nanmean(clus(ir,:));
%     cave=cave+1;
for ir=1:size(clus,1)
    if any(clus(ir,:)==invalidValue)
        toremove(c)=ir;
        c=c+1
    end
    if smoothN>0
        clus(ir,:)=medfilt1(clus(ir,:),smoothN);
    end
end

if removeProfileWithInvalid
    clus(toremove,:)=[];
end
numProfilesFiltered=size(clus)

%see all the material data for the selection
h=figure
if ~seePlots
    set(h,'Visible','off')
end

cm=loadCmap(color_cluster{DescriptorIndex},128);
colormap(cm)
%colormap('jet')
imagesc(clus)
title([selCtag tagDesc '-material-profile']);
ylabel('Material points')
xlabel('Time');
set(gca,'XTick',time_axis)
set(gca,'XTickLabel',time_axis) 
if savePlots
    saveas(h,[StatsPath dataset tagfolder filesep 'clusterRaw' tagDesc '-' selCtag '.png'],'png')
end

%%%%%%%%%%%%%%%%%% TREE %%%%%%%%%%%%%%%%%%%%%
% if DescriptorsSet>3
%     D= pdist(clus(:,2:end), distanceT);
% else
    D= pdist(clus, distanceT);   
%end
size(D)
%Ds=squareform(D);
L= linkage(D, method);
%cophenet
co=cophenet(L,D)
%Dendogram
I=inconsistent(L);
set(0,'RecursionLimit',2000)
h= figure
if ~seePlots
    set(h,'Visible','off')
end
dendrogram(L)
set(0,'RecursionLimit',2000)
file_dendo=[StatsPath dataset tagfolder filesep tagDesc seltag '_dendogram']
saveas(h,[file_dendo '.png'],'png')

size(I)

%%%%%%%%%%%%%%%%%%%%% CLASSES %%%%%%%%%%%%%%%%%%%%%%%
%C=cluster(L,'cutoff',1);
%mxCluster=6
C=cluster(L,'maxclust',mxCluster);
classes=sort(unique(C))
size(C)

if saveClusterData
    'Save intermediate'
    %save(
end

%export2movit=0
if cluster2movit
    fid=fopen(nameSel_clus,'w');
    fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s\n','id_center','selection','timestep','x','y','z','id_mother','validation');
end

%%%% Use C to class ClusterData
%time / spatial_id / descriptorValue
for cl=1:length(classes)
    Class=ClusterData(C==cl,:,:);
    'Class size'
    cl
    sclass=size(Class,1)
    if DescriptorsSet>3
        CImage=clus(C==cl,2:end);
    else
        CImage=clus(C==cl,:);
    end
    
    h=figure
    if ~seePlots
        set(h,'Visible','off')
    end
    colormap(cm)
    imagesc(CImage);
    title([selCtag tagDesc '-class-' num2str(cl)]);
    ylabel('Material points')
    xlabel('Time');
    set(gca,'XTick',[1:length(time_axis):size(CImage,2)])  
    set(gca,'XTickLabel',time_axis) 
    if savePlots
        saveas(h,[StatsPath dataset tagfolder filesep 'cluster' tagDesc '-' selCtag '-class-' num2str(cl) '.png'],'png')
    end
    
    %cluster2movit=0
    if cluster2movit
        %we relabel the spatial ids!
        class_sel=zeros(sclass*size(ClusterData,2),9);%spatio-temporal data again
        class_sel(:)=NaN;
        cd=1
        for ci=1:sclass
            for ti=1:size(ClusterData,2)
                %to come back to the selection of cells we need the
                %spatial/cell id
                spat_id=Class(ci,ti,2);
                data=selec(selec(:,1)==spat_id,:);
                %Replace the original selection number by the clusterized
                %class cl
                if size(data,1)==1
                    %data(1,2)=cl;
                    %nameSel_clus
                    %class_sel(cd,:)=data;
                    %cd=cd+1;
                    fprintf(fid,'%i;%i;%i;%i;%i;%i;%i;%i\n',data(1,1),cl,data(1,3),data(1,4),data(1,5),data(1,6),data(1,7),data(1,8));
                elseif size(data,1)==0
                    'Linked outside the selection'
                    %the lagrangian tracking linked to a cell outside the
                    %selection
                    if add2selection
                        %class_sel(cd,:)=[spat_id cl ti 0 0 0 0 0 0]%The rest is NaNs (hopefully it works)
                        fprintf(fid,'%i;%i;%i;%i;%i;%i;%i;%i\n',spat_id,20,ti,0,0,0,0,0);
                        %cd=cd+1
                    end
                elseif size(data,1)>1
                    'Error'
                    pause
                end
                %pause
            end
        end
    end
end

if cluster2movit
    fclose(fid);
end





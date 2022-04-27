% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%%%%%%%%%%%%%%%%% SPATIO-TEMPORAL CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%
addpath('colormaps')

%%%% Load coherent trajectories
%track_file=[desc 'tracking-lag-' num2str(tinilag) '-' num2str(tfinal) '.csv']
track_file=[desc 'tracking-' num2str(tinilag) '-' num2str(tfinal) xtag '.csv']
Track_lag=dlmread(track_file);
lag_steps=[max(t_start,tinilag):min(t_limit,tfinal)];%It starts in the points we want
mat_ids=[min(Track_lag(:,1)):max(Track_lag(:,1))];%The ids
mat_ids=[unique(Track_lag(:,1))];
num_ids=length(mat_ids)

if old
    t_col=3;
else
    t_col=6;
end
%Careful, here we consider consistent ids! Otherwise we should select by
%time
t_ref=tinilag
if only_selection
    'Restrict Data to the Selections'
    %do the clustering only for the selections u want
    %if a temporal reference is given, the ids are selected in there
    if t_ref>-1
        selec_t=selec(selec(:,3)==t_ref,:);
        Track_lag_t=Track_lag(Track_lag(:,t_col)==t_ref,:);
        Desc_t=Desc(Desc(:,3)==t_ref,:);
    else
        Track_lag_t=Track_lag;
    end
    
    %see the spatial ids for the selectionss
    selec_i= arrayfun(@(x) find(selec_t(:,2)==x), ss, 'UniformOutput', false);
    size(selec_i)
    selec_ii=[]
    for i=1:length(selec_i)
        selec_ii=vertcat(selec_ii,selec_i{i});
    end
    
    %look in the tracking to get the corrsponing material_ids
    spatial_ids=selec_t(selec_ii,1);
    length(spatial_ids);
    spatial_ids=unique(spatial_ids)
    length(spatial_ids)

    %take the rows from Desc of the material ids within the selections
    selec2= arrayfun(@(x) find(Track_lag_t(:,2)==x), spatial_ids, 'UniformOutput', false);
    size(selec2)

    selec2= cell2mat(selec2);
    sTrack_lag=Track_lag_t(selec2,:);
else
    sTrack_lag=Track_lag;
end

material_ids=unique(sTrack_lag(:,1));
num_ids=length(material_ids)

%One dimension for the cluster
%mat_ids * steps * {time/spatial_id/desc}
ClusterData=zeros(num_ids,(lag_steps(end)-lag_steps(1)+1),2+1);
ClusterData(:)=NaN;%we want only complete trajectories

%Crop in time
t_start=tinilag
Track_lag=Track_lag((Track_lag(:,t_col)>=t_start) & (Track_lag(:,t_col)<=t_limit),:);

%Go through all the material data to create the Cluster Datas
toremove=[];
ce=0;
for ii=1:num_ids
    'Material id profile'
    m_id=material_ids(ii)
    %take the trajectory->spatio-temporal data
    traj=Track_lag(Track_lag(:,1)==m_id,[2 t_col]);

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
            data=Desc(selecids,[3 1 index(DescriptorIndex)]);
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
            data=Desc(selecids,[3 1 index(DescriptorIndex)]);
            ClusterData(ii,:,:)=data;
        end
    end
end
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
figure;
c=ClusterData(:,:,3);
c(c>prctile(c(:),98))=prctile(c(:),98);
c(c<prctile(c(:),2))=prctile(c(:),2);
surf(c)
colormap('gray')
axis off

for tcurrent=tinilag:tfinal
    Im=ClusterData(:,tcurrent-tinilag+1,3);

    mx=prctile(Im, pmax(DescriptorIndex));
    mn=prctile(Im, pmin(DescriptorIndex));
    if removeOutliars
        Im=Im(Im>=mn);
        Im=Im(Im<=mx);
    else
        Im(Im>mx)=mx;
        Im(Im<mn)=mn;
    end
    samplesInside=length(Im);
    
    metric(tcurrent+1,1,si)=mean(Im(:));
    metric(tcurrent+1,2,si)=median(Im(:));
    metric(tcurrent+1,3,si)=std(Im(:));
    metric(tcurrent+1,4,si)=mx;%max and min
    metric(tcurrent+1,5,si)=mn;
    metric(tcurrent+1,6,si)=prctile(Im, 99);
    metric(tcurrent+1,7,si)=prctile(Im, 1);
    metric(tcurrent+1,8,si)=prctile(Im, 75);
    metric(tcurrent+1,9,si)=prctile(Im, 25);
    metric(tcurrent+1,10,si)=size(ClusterData,1);
    metric(tcurrent+1,11,si)=size(ClusterData,1);
    metric(tcurrent+1,12,si)=size(ClusterData,1);
    metric(tcurrent+1,13,si)=samplesInside;%data after outliars
    metric(tcurrent+1,14,si)=ss;%the actual selection number
    metric(tcurrent+1,15,si)=tcurrent;%the actual selection number
end


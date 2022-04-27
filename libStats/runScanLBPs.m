% /* --------------------------------------------------------------------------------------
%  * File:    runScanLBP.m
%  * Date:    01/06/2015
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2013-2017, David Pastor Escuredo
%  with Biomedical Image Technology, UPM (BIT-UPM)
%  with BioEmergences, CNRS
%  with LifeD lab
%  All rights reserved.

%%%%%%%%%%%%%%%%% SPATIO-TEMPORAL CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%
addpath('colormaps')

if clus_method
    tag_modes=['-hierarchical' '-' tagDis]
else
    tag_modes=['-kmeans' '-' tagDis]
end

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
    selec_i= arrayfun(@(x) find(selec_t(:,2)==x), clusterSelecBasis, 'UniformOutput', false);
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
ClusterData=zeros(num_ids,(lag_steps(end)-lag_steps(1)+1),2+length(DescriptorsIndexes));
ClusterData(:)=NaN;%we want only complete trajectories

%Crop in time
Track_lag=Track_lag((Track_lag(:,t_col)>=t_start) & (Track_lag(:,t_col)<=t_limit),:);

if redo
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
            data=Desc(selecids,[3 1 index(DescriptorsIndexes)]);
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
            data=Desc(selecids,[3 1 index(DescriptorsIndexes)]);
            ClusterData(ii,:,:)=data;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time constrains %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
 ClusterData(toremove,:,:)=[];
    size(ClusterData)
    save([desc 'ClusterData' seltagBasis],'ClusterData')
else
    load([desc 'ClusterData' seltagBasis])
    size(ClusterData)
end

min_t=min(ClusterData(1,:,1));
max_t=max(ClusterData(1,:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CREATING BASIS FOR EACH DESCRIPTOR %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mxp=98
mnp=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CREATING BASIS FOR EACH DESCRIPTOR %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs=1
for DescriptorIndex=DescriptorsIndexes
    
    clus=squeeze(ClusterData(:,:,cs+2));%just the data but we keep ids
    cs=cs+1
    tagDesc=['-' tagsCmaps{DescriptorIndex} tagg]
    nameDesc=tags{DescriptorIndex};
    nameDescLeg=tagsShort{DescriptorIndex};
    if redo
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering profiles together %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Remove invalid outliars in a material perspective that will affect the clustering
        %Filtering to smooth if necessary
        toremove=[];
        c=1;
        numOfMaterialProfiles=size(clus)
        %general statistical cropping
        mx=prctile(clus(:),mxp)
        mn=prctile(clus(:),mnp)
        clus(clus>mx)=mx;
        clus(clus<mn)=mn;
        %     ave(cave)=nanmean(clus(ir,:));
        %     cave=cave+1;
        clus_old=clus;
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
        
        %%%%%%%%%%%%%% Plotting Material Profile
        %after smoothing
        %see all the material data for the selection
        h=figure
        if ~seePlots
            set(h,'Visible','off')
        end
        %Set the colormap by default
        cm=loadCmap(tagsCmaps{DescriptorIndex});
        colormap(cm)
        % climits=desc_limits;
        %Set the limits of the color!!
        if isSymmetric(DescriptorIndex)
            if DescriptorsSet>3
                clus2=log(clus_old);
                nameDescP=['log(' nameDesc ')'];
            else
                nameDescP=nameDesc ;
                clus2=clus_old;
            end
            maxrange=max(abs([prctile(clus2(:),30) prctile(clus2(:),70)]));
            rangeD=[-maxrange maxrange];
            if use_defaultRange
                rangeD=climits(DescriptorIndex,:)
                rangeD=[-0.01 0.01]
            end
            imagesc(clus2,rangeD);
            
        else
            rangeD=[0 prctile(clus_old(:),90)];
            if use_defaultRange
                rangeD=climits(DescriptorIndex,:)
            end
            imagesc(clus_old,rangeD);
            nameDescP=nameDesc;
        end
        %Plot
        %Legend and texts
        if textSimple==0
            ylabel('Material points','FontSize',fsize);
            title([selCtag tagDesc '-material-profile'],'FontSize',fsize+4);
        elseif textSimple==1
            ylabel(['Trajectories: ' num2str(size(clus,1))],'FontSize',fsize);
            title([nameDesc],'FontSize',fsize+4);
            set(gca,'YTickLabel',[]);
            %elseif textSimple==1
            %    title([tagDesc],'FontSize',fsize+4);
        elseif textSimple==2
            set(gca,'YTickLabel',[]);
            set(gca,'YTick',[]);
        end
        
        %%%%%%%%%%%%%% TIME AXIS %%%%%%%%%%%%%%%%%%%%%%
        dis_time=[min_t:step_time:max_t];
        time_axis=dis_time;
        t_label='Timestep'
        if phys_time
            time_axis=dis_time*dt_min+t0_min;
            time_axis=time_axis/60;
            t_label='Time (hours)'
        end
        
        %set(gca,'XTick',[1:step_time:size(ClusterData,2)])
        % set(gca,'XTickLabel',time_axis)
        if textSimple<2
            xlabel(t_label,'FontSize',fsize)
            set(gca,'XTick',[1:step_time:size(ClusterData,2)])
            set(gca,'XTickLabel',time_axis)
        else
            set(gca,'XTickLabel',[])
            set(gca,'XTick',[])
        end
        if plotHPF
            hpf_time=ceil(([hpf8 hpf9 hpf10 hpf11 hpf12 hpf13]-t0_min)/dt_min)-tinilag+1;
            set(gca,'XTick',hpf_time)
            set(gca,'XTickLabel',{'8','9','10','11','12','13'},'FontSize',fsize)
            if textSimple<2
                xlabel('Time (hpf)','FontSize',fsize)
            end
        end
        colorbar('FontSize',fsize)
        %Big Screen
        screen_size = get(0, 'ScreenSize');
        set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        if savePlots
            saveas(h,[StatsPath dataset tagfolder filesep 'materialProfile' tagDesc '-' selCtag '.png'],'png')
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save([desc 'profile-' tagDesc '-' selCtag], 'clus');    
    size(clus)
    m_stats=[]
    c=1
    for i=1:ww:size(clus,2)-wd
       clusC=clus(:,i:i+wd);
       x=double([mean(clusC(:)) std(clusC(:)) entropy(clusC(:))])
       m_stats(c,:)=x;
       c=c+1
    end
    figure
    plot(m_stats)
    pause
end

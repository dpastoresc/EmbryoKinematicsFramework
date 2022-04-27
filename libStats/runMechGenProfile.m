% /* --------------------------------------------------------------------------------------
%  * File:    runMechGenProfile.m
%  * Date:    01/06/2015
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2013-2018, David Pastor Escuredo
%  with Biomedical Image Technology, UPM (BIT-UPM)
%  with BioEmergences, CNRS
%  All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%% SPATIO-TEMPORAL CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%
addpath('colormaps')
addpath('../../matlab_code/io/')
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
    size(spatial_ids);
    spatial_ids_unique=spatial_ids(:,1);
    spatial_ids_unique=unique(spatial_ids_unique);
    length(spatial_ids_unique);
    %take the rows from Desc of the material ids within the selections
    selec2= arrayfun(@(x) find(Track_lag_t(:,2)==x), spatial_ids_unique, 'UniformOutput', false);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% READING GENETIC PROFILE %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GenData=double(zeros(num_ids,(lag_steps(end)-lag_steps(1)+1),2+3));%max and mean for every gene sample
GenData(:)=NaN;%we want only complete trajectories
mm=3
if redo
    for ti=tinilag:t_max_cluster
        tr=ti-tinilag+1
        stamp_i=makeTimeStamp(ti);
        firstimage=[folder 'VTK/' d stamp_i chg format]
        [Im sin sp]=loadImage(firstimage);
        [Im iinf]=ReadData3D(firstimage);
        traj_t=Desc(Desc(:,3)==ti,:);
        for ii=1:num_ids
            m_id=material_ids(ii);
            p_t=traj_t(traj_t(:,2)==m_id,:);
            
            if size(p_t,1)>0
            cid=p_t(1);
            ll=p_t(4:6)./sp;
            ll=uint16(ll);
            
            xm1=ll(1)-mm;
            xm1(xm1<1)=1;
            xm2=ll(1)+mm;
            xm2(xm2>sin(1))=sin(1);
            
            ym1=ll(2)-mm;
            ym1(ym1<1)=1;
            ym2=ll(2)+mm;
            ym2(ym2>sin(2))=sin(2);
            
            zm1=ll(3)-mm;
            zm1(zm1<1)=1;
            zm2=ll(3)+mm;
            zm2(zm2>sin(3))=sin(3);
            
            v=Im(xm1:xm2,ym1:ym2,zm1:zm2);
            v1=double(median(v(:)));
            v2=double(prctile(v(:),98));
            data=[ti cid v1 v2 m_id];
            GenData(ii,tr,:)=data;
            else
            GenData(ii,tr,:)=[ti NaN NaN NaN m_id];  
            end
        end
    end
    
    %Go through all the material data to create the Cluster Data
    toremove=[];
    ce=0;
    for ii=1:num_ids
        
        'Material id profile'
        m_id=material_ids(ii);
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
    GenData(toremove,:,:)=[];
    if t_max_cluster>-1
        itcluster=find(ClusterData(1,:,1)==t_max_cluster)
    else
        itcluster=100;
    end
    
    ClusterData=ClusterData(:,1:itcluster,:);
    GenData=GenData(:,1:itcluster,:);
    
    ce
    num_ids
    size(ClusterData)
    min_t=min(ClusterData(1,:,1));
    max_t=max(ClusterData(1,:,1));
    
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
    GenData(toremove,:,:)=[];
    cs=1;
    size(ClusterData)
    size(GenData)
    if cN
        save([desc 'GenDataNuc' seltagBasis],'GenData')
        save([desc 'ClusterDataNuc' seltagBasis],'ClusterData')
    elseif cN==0
        save([desc 'GenData' seltagBasis],'GenData')
        save([desc 'ClusterData' seltagBasis],'ClusterData')
    end
else
    if cN==1
        load([desc 'GenDataNuc' seltagBasis])
        load([desc 'ClusterDataNuc' seltagBasis])
    elseif cN==2
        load([desc 'GenDataNuc' seltagBasis])
        load([desc 'ClusterDataNuc' seltagBasis])
        GC=GenData;
        load([desc 'GenData' seltagBasis]);
        load([desc 'ClusterData' seltagBasis]);
        %GenData=GenData./GC;
    elseif cN==3
        load([desc 'GenDataNuc' seltagBasis])
        load([desc 'ClusterDataNuc' seltagBasis])
        GC=GenData;
        load([desc 'GenData' seltagBasis]);
        load([desc 'ClusterData' seltagBasis]);    
    else
        load([desc 'GenData' seltagBasis])
        load([desc 'ClusterData' seltagBasis])  
    end
    size(ClusterData)
    size(GenData)
end
cs=1
if genclus
    ClusterData=GenData;
    modes_colors={[ 227 237 0; 194 199 10; 200 190 20;],[227 237 0; 194 199 10; 200 190 20;],[227 237 0; 194 199 10; 200 190 20;],[227 237 0; 194 199 10; 200 190 20;]}
    cs=2;
    %DescriptorsIndexes=1
end
min_t=min(ClusterData(1,:,1));
max_t=max(ClusterData(1,:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CREATING BASIS FOR EACH DESCRIPTOR %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mxp=98
mnp=2

if cN>1
clusNUC=squeeze(GC(:,:,cs+2));
end
for DescriptorIndex=DescriptorsIndexes
    
    clus=squeeze(ClusterData(:,:,cs+2));%just the data but we keep ids
    cs=cs+1
    if genclus
        if cN==1
            tagDesc='-nuc'
             nameDesc='nuc'
             nameDescLeg='nuc'
        elseif cN==0
            tagDesc='-gsc'
            nameDesc='gsc'
            nameDescLeg='gsc'
        elseif cN==2
            tagDesc='-gscCor'
            nameDesc='gscCor'
            nameDescLeg='gscCor'
        elseif cN==3
            tagDesc='-gscN'
            nameDesc='gscN'
            nameDescLeg='gscN'
        end
        cs=200
        mxp=99
        mnp=0
        smoothN=5
        pgc=100
    else
        tagDesc=['-' tagsCmaps{DescriptorIndex} tagg]
        nameDesc=tags{DescriptorIndex};
        nameDescLeg=tagsShort{DescriptorIndex};
    end
    
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
    if ~genclus
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
               % rangeD=[-0.01 0.01]
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
    else
        cm=loadCmap('gsc');
        colormap(cm)
        rangeD=[0 prctile(clus_old(:),pgc)];
        imagesc(clus_old,rangeD);
        nameDescP=nameDesc;
    end
    
    %%%%% find NAN in gene%%%%
    if genclus
        c=1
        for ir=1:size(clus,1)
            if any(isnan(clus(ir,:)))
                toremove(c)=ir;
                c=c+1
            end
        end
        size(toremove)

        clus(toremove,:)=[]
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
       % saveas(h,[StatsPath dataset tagfolder filesep 'materialProfile' tagDesc '-' selCtag '.png'],'png')
    end
    
    %%%%%%%%%%%%%%%%%% UNSUPERVISD CLUSTERING %%%%%%%%%%%%%%%%%%%%%%%%
    %Using K-means for analysis
    if cN==2
            clus=clus./clusNUC;
    end
    if useDefault
        mxCluster=modes_mech(DescriptorIndex)
    end
    if clus_method==0
        if startCluster
            [C ctr]=kmeans(clus, mxCluster, 'Distance', distanceK, 'replicates', 20, 'start', 'cluster')
        else
            [C ctr]=kmeans(clus, mxCluster, 'Distance', distanceK, 'replicates', 20)
        end
        %size(idx)
        size(ctr)
        cn=isnan(C);
        C(cn)=mxCluster+1;
        classes=sort(unique(C))
    elseif clus_method==1
        D= pdist(clus, distanceT);
        %end
        size(D)
        %Ds=squareform(D);
        L= linkage(D, methodHier);
        %cophenet
        co=cophenet(L,D)
        %Dendogram
        I=inconsistent(L);
        set(0,'RecursionLimit',reclimit)
        
        if extraPlots
            h= figure
            if ~seePlots
                set(h,'Visible','off')
            end
            dendrogram(L)
            file_dendo=[StatsPath dataset tagfolder filesep 'Modes-Tree' tagDesc seltagBasis]
            saveas(h,[file_dendo '.png'],'png')
        end
        size(I)
        
        %%% Get Classes
        %C=cluster(L,'cutoff',1);
        %mxCluster=6
        C=cluster(L,'maxclust',mxCluster);
        cn=isnan(C);
        C(cn)=mxCluster+1;
        classes=sort(unique(C));
        size(C)
    end
    
    %%%%%%%%%%%%%%%%% EXTRACTING MAIN SIGNAL FOR MODES
    %%%% Each class is a signal!
    classes=classes(1:mxCluster);
    lenBasis=length(classes);%NANS OUT
    Basis=zeros(length(classes),size(clus,2));
    DataB=zeros(3,lenBasis);
    V=zeros(size(clus,1), lenBasis);
    V(:)=NaN;
    if cN==3
            clus=clus./clusNUC;
    end
    %%%%%%%%%%%%%%%%%%%% Get Main Signal of the MODE %%%%%%%%%%%%%%%%%%
    for cl=1:lenBasis
       
        CImage=clus(C==cl,:);
        sclass=size(CImage,1)
        signal=mean(CImage);%MODE AS THE MEAN!!!!
        DisBase=zeros(size(clus,1),1);
        DisBase(:)=NaN;
        for sss=1:sclass
            if distanceType==1
                DisBase(sss)=pdist2( CImage(sss,:), signal,'euclidean');
            elseif distanceType==2
                DisBase(sss)=pdist2( CImage(sss,:), signal,'cosine');
            elseif distanceType==0
                'Correlation'
                DisBase(sss)=pdist2( CImage(sss,:), signal,'correlation');
            end
        end
        %V(:,cl)=DisBase;
        Basis(cl,:)=signal;
        %DataB(1,cl)=nanmean(DisBase);
        %DataB(2,cl)=nanstd(DisBase);
        %DataB(3,cl)=sclass;
        %Ordered=vertcat(Ordered, CImage);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Reorder profile by clustering %%%%%%%%%%%%%%%%%%%
    mbasis=median(Basis,2)
    [mbasis msort]=sort(mbasis,'descend')
    BasisOrdered=zeros(size(Basis));
    Ordered=[];
    for cl=1:lenBasis
        co=msort(cl);
        CImage=clus(C==co,:);
        BasisOrdered(cl,:)=Basis(co,:);
        %%%%%%%%%%%%%%%%% Reorder profile by clustering %%%%%%%%%%%%%%%%%%%
        Ordered=vertcat(Ordered, CImage);
    end
    Basis=BasisOrdered;%replace for an ordered version
    %%%% SAVE MODES BASIS
    save([desc 'Basis' tagDesc seltagBasis tag_modes '-modes' num2str(mxCluster)],'Basis')
    %load([desc 'Basis' tagDesc seltagBasis tag_modes '-modes' num2str(mxCluster)])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Plotting Material Profile after reordering %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h=figure
    if ~seePlots
        set(h,'Visible','off')
    end
    if ~genclus
    %Set the colormap by default
    cm=loadCmap(tagsCmaps{DescriptorIndex});
    colormap(cm)
    %Set the limits of the color!!
    if isSymmetric(DescriptorIndex)
        clus2=log(Ordered);
        maxrange=max(abs([prctile(clus2(:),30) prctile(clus2(:),70)]));
        rangeD=[-maxrange maxrange];
        if use_defaultRange
            rangeD=climits(DescriptorIndex,:)
        end
        imagesc(clus2,rangeD);
        nameDescP=['log(' nameDesc ')'];
    else
        rangeD=[0 prctile(clus_old(:),90)];
        if use_defaultRange
            rangeD=climits(DescriptorIndex,:)
        end
            imagesc(Ordered,rangeD);
            nameDescP=nameDesc;
    end
    else
        cm=loadCmap('gsc');
        colormap(cm)
        rangeD=[0 prctile(clus_old(:),pgc)];
        imagesc(Ordered,rangeD);
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
        saveas(h,[StatsPath dataset tagfolder filesep 'materialProfileOrdered' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
    end
   
    figure
    if clus_method
        [s,h]=silhouette(clus, C, distanceT);
    else
        [s,h]=silhouette(clus, C, distanceK);
    end
    
    %set(gca,'Color',[.8 .8 1]);
    set(gca,'FontSize',fsize)
    sils=allchild(gca);
    %for sll=1:length(sils)
    set(sils(6),'FaceColor',[.8 .8 1]);
    %end
    %set(gca,'YTickLabel',{'1','2','3','4','5'},'FontSize',fsize)
    xlabel('Silhouette value','FontSize',fsize)
    ylabel('Mode number','FontSize',fsize)
    screen_size = get(0, 'ScreenSize');
    set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    
    if savePlots
        saveas(h,[StatsPath dataset tagfolder filesep 'modesSil' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT MODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h=figure
    hold on
    if ~seePlots
        set(h,'Visible','off')
    end
    %PLOT
    if useDefault
        for icc=1:size(Basis,1)
            if ~genclus
                colorrgb=modes_colors{DescriptorIndex}(icc,:)/255
            else
                if cN~=1
                modes_colors={[255 21 10;255 194 79;65 255 255;255 63 255]}
                modes_colors={[21 255 10;65 255 255;255 0 255;   255 63 255]}
                colorrgb=modes_colors{1}(icc,:)/255     
                else
                modes_colors={[0.9 0.9 0.9;0.7 0.7 0.7;0.5 0.5 0.5;0.3 0.3 0.3]}
                colorrgb=modes_colors{1}(icc,:);
                end
            end
            plot(Basis(icc,:)','LineWidth',widthMode,'Color',colorrgb);
            %pause
        end
    else
        %modes_colors={[255 21 10; 255 194 79; 255 63 255]}
        %colorrgb=modes_colors{1}(icc,:)/255 
        plot(Basis(icc,:)','LineWidth',widthMode)%,'Color',colorrgb);
    end
    if gridOn
        grid on
    end
    %Title
    if textSimple<1
        if clus_method==0
            title([selCtag tagDesc '- EigenProcesses (K-means ' distanceK ' )'],'FontSize',fsize+4);
        else
            title([selCtag tagDesc '- EigenProcesses (Hierarchical ' distanceT ' )','FontSize',fsize+4]);
        end
    end
    
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
    if ~genclus
    if use_defaultRange
        set(gca,'YLim', desc_limits(DescriptorIndex,:));
    end
    end
    %Legend
    for ib=1:size(Basis,1)
        mm_lg{ib}=[nameDescLeg '_{m' num2str(ib) '}'];
    end
    legend(mm_lg,'Location','NorthOutside','Orientation','horizontal','FontSize',fsize)
    %Save
    screen_size = get(0, 'ScreenSize');
    set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    hold off
    if savePlots
        if useDefault
          %  saveas(h,[StatsPath dataset tagfolder filesep 'modesProfileSelected' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
            if cN==1
                saveas(h,[StatsPath dataset tagfolder filesep 'modesProfileSelectedNUC' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
            elseif cN==2
                saveas(h,[StatsPath dataset tagfolder filesep 'modesProfileSelectedCorr' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')   
            elseif cN==3
                saveas(h,[StatsPath dataset tagfolder filesep 'modesProfileSelectedCorr2' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')   
            else 
                'saving'
                saveas(h,[StatsPath dataset tagfolder filesep 'modesProfileSelected' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
            end
        else
            if cN==1
                saveas(h,[StatsPath dataset tagfolder filesep 'modesProfileNUC' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
            elseif cN==2
                saveas(h,[StatsPath dataset tagfolder filesep 'modesProfileCorr' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')   
            elseif cN==3
                saveas(h,[StatsPath dataset tagfolder filesep 'modesProfileCorr2' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')   
            else
                saveas(h,[StatsPath dataset tagfolder filesep 'modesProfile' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROFILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gene is the classifier
    
    if genProfiling
        load([desc 'ClusterData' seltagBasis])
        mindp=[0 0 0]
        maxdp=[0 0 0]
        for cl=1:mxCluster
            co=msort(cl);
            Class=ClusterData(C==co,:,:);
            ccc=1
            for di=[3 4 6]
                descript=Class(:,:,di)
                dp=mean(descript)
                mindp1=min(dp);
                maxdp1=max(dp);
                mindp(ccc)=min(mindp(ccc),mindp1);
                maxdp(ccc)=max(maxdp(ccc),maxdp1);
                ccc=ccc+1
            end
        end
        for cl=1:mxCluster
            co=msort(cl);
            h=figure
            hold on
            bp=Basis(cl,:);
            minbp=min(Basis(:));
            maxbp=max(Basis(:));
            %bp=der1D(b)
            if ac
                bp=ac1D(bp); 
            end
            if normalizeG
                bp=(bp-minbp)./(maxbp-minbp);
            end
            modes_colors={[255 21 10;255 194 79;65 255 255;   255 63 255]}
            modes_colors={[21 255 10;65 255 255;255 0 255;   255 63 255]}
            crgb=modes_colors{1}(cl,:)/255
            plot(bp,'LineWidth',widthMode,'Color',crgb)
            Class=ClusterData(C==co,:,:);
            
            crgb=[0.9 0.9 0.9; 0.7 0.7 0.7; 0.5 0.5 0.5]
            ccol=1
            for di=[3 4 6]
                crgb1=crgb(ccol,:);
                
                descript=Class(:,:,di)
                dp=mean(descript)
                if der
                    dp=der1D(dp)
                end
                
                dp=((dp-mindp(ccol))./(maxdp(ccol)-mindp(ccol)))
                plot(dp,'LineWidth',widthMode,'Color',crgb1)
                ccol=ccol+1;
            end
            if textSimple<1
                title([selCtag tagDesc '-genevsdeformation','FontSize',fsize+4]);
            end
            
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
            mm_lg={'gsc','\Delta V','\Delta \gamma _1','\Delta \alpha'}
            legend(mm_lg,'Location','NorthOutside','Orientation','horizontal','FontSize',fsize)
            %Save
            
            set(gca,'YLim',[0 1])
 
            hold off
            screen_size = get(0, 'ScreenSize');
            set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            if savePlots
                saveas(h,[StatsPath dataset tagfolder filesep 'segmentation-profile' tagDesc '-' selCtag '-seg' num2str(cl) '.png'],'png')
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXPORT MODES CLUSTERS %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sss=0
    add2selection=0
    if cluster2movit
        nameSel_modes=[datapath dataset tagfolder '_t' filesep 'modesClusters' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) xtag '-' num2str(tinilag)  '-' num2str(t_max_cluster) '.csv']
        fid=fopen(nameSel_modes,'w');
        fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s\n','id_center','selection','timestep','x','y','z','id_mother','validation');
        
        for cl=1:mxCluster
            co=msort(cl);
            Class=ClusterData(C==co,:,:);
            'Class size'
            cl
            size(Class)
            sclass=size(Class,1)
            sss=sss+sclass
            %we relabel the spatial ids!
            class_sel=zeros(sclass*size(ClusterData,2),9);%spatio-temporal data again
            class_sel(:)=NaN;
           
            for ci=1:sclass
                for ti=1:1%size(ClusterData,2)
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
                        %'Linked outside the selection'
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
        fclose(fid);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXPORT METRIC MAP %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exportBasisMap
        
        %Calculate distances to the signal
        DisMaterial=zeros(size(clus,1),size(Basis,1));
        for lb=1:size(Basis,1)
            for lc=1:size(clus,1)
                distc=sum((clus(lc,:)-Basis(lb,:)).^2).^0.5
                DisMaterial(lc,lb)=distc;
            end
        end
        
        maxBasisDis=max(DisMaterial);
        minBasisDis=min(DisMaterial);
        
        ll=size(ClusterData,2)*size(DisMaterial,1)
        
        DescriptorsM=zeros(ll,12)
        lcount=1
        for ti=1:size(ClusterData,2)
            Cluster_t=ClusterData(:,ti,:);
            size(Cluster_t)
            for ld=1:size(DisMaterial,1)
                %DisMaterial(ld,:)
                cellid=Cluster_t(ld,2);
                DescEntry=Desc(Desc(:,1)==cellid, 1:6);%info
                DescEntry=horzcat(DescEntry(1,:), DisMaterial(ld,:));
                DescriptorsM(lcount,:)=DescEntry;
                lcount=lcount+1;
            end
        end
        MechanomeDescriptors2movit;
    end
end

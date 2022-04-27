% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%%%%%%%%%%%%%%%%% SPATIO-TEMPORAL CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%
'Projecting tissue'
addpath('colormaps')
if clus_method
    tag_modes=['-hierarchical' '-' tagDis]
else
    tag_modes=['-kmeans' '-' tagDis]
end

if segmentField
    nameSel_mechanome=[datapath dataset tagfolder '_t' filesep dataset '_mechseg_t' num2str(tinilags(dN)) '-tf' num2str(t_max_cluster)...
        seltag '-' num2str(length(DescriptorsIndexes)) '-modes' num2str(segmentationClasses) tag_modes...
        '-' distanceClass '-' method '-' num2str(sameRangeDis) '-' num2str(cropDisRange) '-' num2str(binarizeProjection) xtag '.csv']
else
    cluster2movit=0
    nameSel_mechanome=''%['mechanome_t' num2str(tinilags(1)) seltag 'modes_seg3.csv']
end

if old
    t_col=3;
else
    t_col=6;
end

%%%% Load coherent trajectories
%track_file=[desc 'tracking-lag-' num2str(tinilag) '-' num2str(tfinal) '.csv']
track_file=[desc 'tracking-' num2str(tinilag) '-' num2str(tfinal) xtag '.csv']
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
        Track_lag_t=Track_lag(Track_lag(:,t_col)==t_ref,:);
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
    spatial_ids=selec_t(selec_ii,[1 2]);
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
'num ids'

%One dimension for the cluster
%mat_ids * steps * {time/spatial_id/desc}
ClusterData=zeros(num_ids,(lag_steps(end)-lag_steps(1)+1),2+length(DescriptorsIndexes));
ClusterData(:)=NaN;%we want only complete trajectories

%Crop in time
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
time_axis=[min_t:step_time:max_t];
%Remove uncomplete trajectories (the time has been already cropped)
toremove=[];
c=1;
numOfOriginalIds=size(ClusterData,1);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate distances to Modes %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DisMaterialAll=[];
cs=1
css=1
m_xlabel=[]
m_color=[]
label_count=1
violin_color=[0 0 0; 0 0 0; 0 0 0; 0 0 0; 0.8 0.4 0.4; 0.4 0.8 0.8; 0.4 0.4 0.8; 0.4 0.8 0.4];

for DescriptorIndex=DescriptorsIndexes
    
    'Processing' 
    DescriptorIndex
    
    clus=squeeze(ClusterData(:,:,cs+2));%just the data but we keep ids
    tagDesc=['-' tagsCmaps{DescriptorIndex} tagg]
    nameDesc=[tags{DescriptorIndex}]
    nameDescLeg=[tagsShort{DescriptorIndex}]
    mxCluster=modes_mech(DescriptorIndex)
    
    %%%%%%%%%%%%%%% Load selected basis
    tag_clus=[tag_modes '-modes' num2str(mxCluster)]
    load([desc 'Basis' tagDesc seltag tag_clus])
    
    %%%%%%%%%%%%%%%% Filter material profiles %%%%%%%%%%%%%%%%%
    %Remove invalid outliars in a material perspective that will affect the clustering
    %Filtering to smooth if necessary
    %ave=[]
    %cave=1
    toremove=[];
    c=1;
    numOfMaterialProfiles=size(clus)
    
    %liitation in number of profiles.
    if limitP>0
        clus=clus(1:limitP,:)
    end
    
    %general statistical cropping
    mx=prctile(clus(:),pmax(DescriptorIndex))
    mn=prctile(clus(:),pmin(DescriptorIndex))
    clus(clus>mx)=mx;
    clus(clus<mn)=mn;
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
    clus_old=clus;
    %Set the limits of the color!!
    if isSymmetric(DescriptorIndex)
        clus2=log(clus_old);
        maxrange=max(abs([prctile(clus2(:),30) prctile(clus2(:),70)]));
        rangeD=[-maxrange maxrange];
        if use_defaultRange
            rangeD=climits(DescriptorIndex,:)
        end
        imagesc(clus2,rangeD);
        nameDesc2=['log(' nameDesc ')'];
    else
        nameDesc2=nameDesc
        rangeD=[0 prctile(clus_old(:),90)];
        if use_defaultRange
            rangeD=climits(DescriptorIndex,:)
        end
        imagesc(clus_old,rangeD);
    end
    %Plot
    %Legend and texts
    if textSimple==0
        ylabel('Material points','FontSize',fsize);taf
        title([selCtag tagDesc '-material-profile'],'FontSize',fsize+4);
    elseif textSimple==1
        ylabel(['Trajectories: ' num2str(size(clus,1))],'FontSize',fsize);
        title([nameDesc2 ' ' selCtag],'FontSize',fsize+4);
        set(gca,'YTickLabel',[]);
        %elseif textSimple==1
        %    title([tagDesc],'FontSize',fsize+4);
    elseif textSimple==2
        set(gca,'YTickLabel',[]);
        set(gca,'YTick',[]);
    end
    
    %%% TIME AXIS
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
    colorbar('FontSize',fsize);
    %Big Screen
    screen_size = get(0, 'ScreenSize');
    set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    if savePlots
        saveas(h,[StatsPath dataset tagfolder filesep 'materialProfile' tagDesc '-' selCtag '.png'],'png')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% DOMAIN PROJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h=figure
    hold on
    if ~seePlots
        set(h,'Visible','off')
    end
    %Plot modes
    for icc=1:size(Basis,1)
        colorrgb=modes_colors{DescriptorIndex}(icc,:)/255
        plot(Basis(icc,:)','LineWidth',widthMode,'Color',colorrgb);        
    end
    if gridOn
        grid on
    end
    %Plot profile of domain
 
    shadedErrorBar([1:size(clus,2)], mean(clus_old),std(clus_old),{'Color',tissue2color(clusterSelec),'LineWidth',widthMode},transp)    

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
    if use_defaultRange
     set(gca,'YLim', desc_limits(DescriptorIndex,:));
    end
    
    %Legend
    mm_lg={}
    for ib=1:size(Basis,1)
        mm_lg{ib}=[nameDescLeg '_{m' num2str(ib) '}'];
        mech_mode_tags{css}=[nameDescLeg '_{m' num2str(ib) '}'];
        css=css+1
    end
    %mm_lg{ib+1}=selCtag
    
    legend(mm_lg,'Location','NorthOutside','Orientation','horizontal','FontSize',fsize)
    %Save
    hold off
    screen_size = get(0, 'ScreenSize');
    set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    
    if savePlots
        saveas(h,[StatsPath dataset tagfolder filesep 'materialProjection' tagDesc '-' selCtag tag_clus '.png'],'png')
    end
    
    %%%%%%%%%%%%%%%%%%% DOMAIN-MECHANOME %%%%%%%%%%%%%%%
    %%%%%%%%%   Calculate distances to each mode    
    if binarizeProjection
        DisMaterialB=zeros(size(clus,1),3);
        DisMaterialC=zeros(size(clus,1),3);
    end
    for lb=1:size(Basis,1)
        DisMaterialO=zeros(size(clus,1),1);
        m_color=vertcat(m_color,modes_colors{DescriptorIndex}(lb,:)/255);
        for lc=1:size(clus,1)
            if distanceType==1
                distc=pdist2( clus(lc,:), Basis(lb,:),'euclidean');
            elseif distanceType==2
                distc=pdist2( clus(lc,:), Basis(lb,:),'cosine');
            elseif distanceType==0
                'Correlation'
                distc=pdist2( clus(lc,:), Basis(lb,:),'correlation');
            end
            DisMaterialO(lc,1)=distc;
        end
        if binarizeProjection
            DisMaterialB(:,lb)=DisMaterialO;
        else
            if sameRangeDis
                if cropDisRange
                    DisMaterialO(DisMaterialO>prctile(DisMaterialO(:),95))=prctile(DisMaterialO(:),95);
                    DisMaterialO(DisMaterialO<prctile(DisMaterialO(:),5))=prctile(DisMaterialO(:),95);
                end
                DisMaterial=max(DisMaterialO(:))-DisMaterialO;
                DisMaterial=DisMaterial./max(DisMaterial(:));
            else
                'no eq'
                DisMaterial=DisMaterialO;
            end
            DisMaterialAll=horzcat(DisMaterialAll, DisMaterial);
        end
        
    end
    if binarizeProjection
        minDesc=min(DisMaterialB,[],2);
        for lb2=1:size(Basis,1)
            check=DisMaterialB(:,lb2)-minDesc;
            DisMaterialC(check==0,lb2)=1;
        end
        DisMaterialAll=horzcat(DisMaterialAll, DisMaterialC);
    end
    %     figure
    %     colormap('Copper')
    %     imagesc(DisMaterial');
    %     title(['Mechanome ' tagDesc])
    %     max(DisMaterial(:))
    %     m_xlabel=horzcat(m_xlabel,['    ' nameDesc]);
    cs=cs+1
end

h=figure
colormap('Copper')
imagesc(DisMaterialAll')
colorbar('FontSize',fsize)
screen_size = get(0, 'ScreenSize');
set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
xlabel(['Trajectories: ' num2str(size(DisMaterialAll,1))],'FontSize',fsize)
% set(gca,'XTickLabel',[],'FontSize',fsize)
% set(gca,'XTick',[])
% rotateXLabels( gca(), 45 )
% set(gca,'YTickLabel',[],'FontSize',fsize)
% set(gca,'YTick',[]])
title('Distances map raw', 'FontSize', fsize)

if segmentField
    %%%%%%%%%%%%%%%%%% TREE %%%%%%%%%%%%%%%%%%%%%
    D= pdist(DisMaterialAll, distanceClass);%the mechanome is semegmented with the distances euclidean vector
    %end
    size(D)
    %Ds=squareform(D);
    L= linkage(D, method);
    %cophenet
    co=cophenet(L,D)
    %Dendogram
    I=inconsistent(L);
    set(0,'RecursionLimit',reclimit)%2000)
    h= figure
    if ~seePlots
        set(h,'Visible','off')
    end
    dendrogram(L)
    set(0,'RecursionLimit',reclimit)%00)
    file_dendo=[StatsPath dataset tagfolder filesep tagDesc seltag '_dendogramMechanome']
    saveas(h,[file_dendo '.png'],'png')
    size(I)
    
    %%%%%%%%%%%%%%%%%%%%% CLASSES %%%%%%%%%%%%%%%%%%%%%%%
    C=cluster(L,'maxclust',segmentationClasses);
    classes=sort(unique(C))
    size(C)
else
    C=ones(size(DisMaterialAll,1),1);
    classes=[1]
end

AllClass=[];
for cl=1:length(classes)
    
    %Statistical distribution of distances
    CImage=DisMaterialAll(C==cl,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% Plot Signature per class %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h=figure
    colormap('Copper')
    imagesc(CImage)
    %colorbar('FontSize',fsize)
    screen_size = get(0, 'ScreenSize');
    set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    %ylabel('Modes','FontSize',fsize)
    %ylabel(['Trajectories: ' num2str(size(CImage,1))],'FontSize',fsize)
    set(gca,'YTickLabel',[],'FontSize',fsize)
    set(gca,'YTick',[])
    set(gca,'XTickLabel',[],'FontSize',fsize)
    set(gca,'XTick',[])
    % set(gca,'XTickLabel',mech_mode_tags,'FontSize',fsize)
    % set(gca,'XTick',[1:size(CImage,2)])
    % rotateXLabels( gca(), 45 )
    
    if savePlots
        saveas(h,[StatsPath dataset tagfolder filesep 'materialMechanome-' selCtag '-' num2str(tinilag) '-modesUsed' num2str(length(DescriptorsIndexes)) '-classes' num2str(length(classes)) '-class' num2str(cl) '.png'],'png')
    end
    
    AllClass=vertcat(AllClass,CImage);
    h=figure
    m_color(m_color>1)=1;
    violin(CImage,'facecolor',m_color,'edgecolor',[],'mc','k', 'medc','k--','bw',0.2)
    %ylabel('Modes of trajectories','FontSize',fsize)
    if textSimple<1
        xlabel(t_label,'FontSize',fsize)
        set(gca,'XTick',[1:size(CImage,2)])
        set(gca,'XTickLabel',[1:size(CImage,2)],'FontSize',fsize)
        if segmentField
            title(['Distance to Mode distribution for cluster ' num2str(cl)],'FontSize', fsize+4);
        else
            title(['Distance to Mode distribution for tissue'],'FontSize', fsize+4);
        end
    else
        set(gca,'XTickLabel',mech_mode_tags,'FontSize',fsize)
        set(gca,'XTick',[1:size(CImage,2)])
        rotateXLabels( gca(), 90 )
        if binarizeProjection
        set(gca,'YTickLabel',{'No', 'Yes'},'FontSize',fsize+2)
        set(gca,'YTick',[0 1])
        end
        %rotateXLabels( gca(), 90 )
    end
    
    taggg=['-t' num2str(tinilag)]
    if segmentField
        taggg=[taggg '-c' num2str(cl)]
    end
    screen_size = get(0, 'ScreenSize');
    set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    grid on
    if savePlots
        saveas(h,[StatsPath dataset tagfolder filesep 'materialProjectionDistribution-' selCtag '-' num2str(tinilag) '-class' num2str(cl) '.png'],'png')
    end
    
    Distrib(cl)=size(CImage,1);
    %Granulometry
    weights=zeros(1,size(CImage,2));
    if granulometry
        for lb=1:size(CImage,2)
            samples=CImage(CImage(:,lb)>0.9,lb);
            size(samples)
            weights(lb)=size(samples,1)/size(CImage,1);
        end
    end
end

if extraPlots
    figure
    bar(weights)
    ylabel('Coefficient - tissue region %','fontsize', 18')
    xlabel('Mode','fontsize', 18)
    figure
    bar(Distrib)
    xlabel('Mechanome clusters','fontsize', 18)
    ylabel('Cells enclosed','fontsize', 18)
    set(gca,'fontsize', 16)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Plot Signature complete %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure
colormap('Copper')
imagesc(AllClass)
%colorbar('FontSize',fsize)
screen_size = get(0, 'ScreenSize');
set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%ylabel('Modes','FontSize',fsize)
%ylabel(['Trajectories: ' num2str(size(DisMaterialAll,1))],'FontSize',fsize)
set(gca,'YTickLabel',[],'FontSize',fsize)
set(gca,'YTick',[])
set(gca,'XTickLabel',mech_mode_tags,'FontSize',fsize)
set(gca,'XTick',[1:size(CImage,2)])
rotateXLabels( gca(), 90 )

if savePlots
    saveas(h,[StatsPath dataset tagfolder filesep 'materialMechanome-' selCtag '-' num2str(tinilag) '-modesUsed' num2str(length(DescriptorsIndexes)) '-classes' num2str(length(classes)) '.png'],'png')
end

%%%%%%%%%%%%% Export 2 movit %%%%%%%%%%%%%%%%%
add2selection=1
if cluster2movit
    fid=fopen(nameSel_mechanome,'w');
    fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s\n','id_center','selection','timestep','x','y','z','id_mother','validation');
end
for cl=1:length(classes)
   
    Class=ClusterData(C==cl,:,:);
    'Class size'
    cl
    sclass=size(Class,1)
    
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
end

if cluster2movit
    fclose(fid);
end

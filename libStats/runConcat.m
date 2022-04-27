% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%%%%%%%%%%%%%%%%% SPATIO-TEMPORAL CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%
addpath('colormaps')

if clus_method
    tag_modes=['-hierarchical' '-' tagDis]
else
    tag_modes=['-kmeans' '-' tagDis]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CREATING BASIS FOR EACH DESCRIPTOR %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs=1


%%%%%%%%%%%%%%%%%% UNSUPERVISD CLUSTERING %%%%%%%%%%%%%%%%%%%%%%%%
%Using K-means for analysis
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
    classes=sort(unique(C))
    size(C)
end

%%%%%%%%%%%%%%%%% EXTRACTING MAIN SIGNAL FOR MODES
%%%% Each class is a signal!
lenBasis=length(classes) ;
Basis=zeros(length(classes),size(clus,2));
DataB=zeros(3,lenBasis);
V=zeros(size(clus,1), lenBasis);
V(:)=NaN;

%%%%%%%%%%%%%%%%%%%% Get Main Signal of the MODE %%%%%%%%%%%%%%%%%%
for cl=1:lenBasis
    
    CImage=clus(C==cl,:);
    sclass=size(CImage,1)
    signal=mean(CImage);%MODE AS THE MEAN!!!!
    signalstd=std(CImage);%MODE AS THE MEAN!!!!
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
    BasisSTD(cl,:)=signalstd;
    %DataB(1,cl)=nanmean(DisBase);
    %DataB(2,cl)=nanstd(DisBase);
    %DataB(3,cl)=sclass;
    %Ordered=vertcat(Ordered, CImage);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Reorder profile by clustering %%%%%%%%%%%%%%%%%%%
mbasis=median(Basis,2)
[mbasis msort]=sort(mbasis,'descend')
BasisOrdered=zeros(size(Basis));
BasisSTDo=zeros(size(Basis));
Ordered=[];
for cl=1:lenBasis
    co=msort(cl);
    CImage=clus(C==co,:);
    BasisOrdered(cl,:)=Basis(co,:);
    %%%%%%%%%%%%%%%%% Reorder profile by clustering %%%%%%%%%%%%%%%%%%%
    Ordered=vertcat(Ordered, CImage);
    BasisSTDo(cl,:)=BasisSTD(co,:);
end
Basis=BasisOrdered;%replace for an ordered version
BasisSTD=BasisSTDo;
%%%% SAVE MODES BASIS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Plotting Material Profile after reordering %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure
if ~seePlots
    set(h,'Visible','off')
end
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
    rangeD=[0 prctile(clus(:),90)];
    if use_defaultRange
        rangeD=climits(DescriptorIndex,:)
    end
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
dis_time=[1:1:120];
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
    set(gca,'XTick',[1:step_time:120])
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
    saveas(h,[StatsPath dataset tagfolder filesep 'materialProfileAllOrdered' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% STATS OF THE MODES (a bit outdated)
%     V=max(V(:))-V;
%     V=V./max(V(:));
%     if extraPlots
%         h=figure
%         violin(V)
%         title(['Modes representativeness distribution ' tagDesc])
%         if savePlots
%             saveas(h,[StatsPath dataset tagfolder filesep 'Modes-distribution' tagDesc '-' selCtag '.png'],'png')
%         end
%         %     figure
%         %     bar(DataB(1,:))
%         %     title('MEAN CLUSTERS->MODES')
%         %     figure
%         %     bar(DataB(2,:))
%         %     title('STD CLUSTERS->MODES')
%         %     figure
%         %     bar(DataB(3,:))
%         %     title('SIZE CLUSTERS->MODES')
%         %     %hold off
%     end

%%% TODO!!
%%% USE pdist + squareform to see the distance between the modes.
%     DisBasis=zeros(lenBasis,lenBasis);
%     for cl=1:lenBasis
%         for ccl=1:lenBasis
%             DisBasis(cl,ccl)=sum( (Basis(cl,:)-Basis(ccl,:)).^2 ).^0.5
%         end
%     end
%     DisBasis=max(DisBasis(:))-DisBasis;
%     DisBasis=DisBasis./max(DisBasis(:));
%     if extraPlots
%         figure
%         imagesc(DisBasis)
%         colormap('Gray')
%     end

%%%% SILHOUETTE
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
    saveas(h,[StatsPath dataset tagfolder filesep 'modesSilAll' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
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
        colorrgb=modes_colors{DescriptorIndex}(icc,:)/255
       % plot(Basis(icc,:)','LineWidth',widthMode,'Color',colorrgb);
        shadedErrorBar(time_axis, Basis(icc,:)',BasisSTD(icc,:)',{'Color',colorrgb},0.8)%,),'LineWidth',tl)
    end
else
    plot(Basis','LineWidth',widthMode);
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
    set(gca,'XTick',[1:step_time:120])
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
set(gca,'XLim',[time_axis(1) time_axis(end)]);
%Legend
for ib=1:size(Basis,1)d
    mm_lg{ib}=[nameDescLeg '_{m' num2str(ib) '}'];
end
%legend(mm_lg,'Location','NorthOutside','Orientation','horizontal','FontSize',fsize)
%Save
screen_size = get(0, 'ScreenSize');
set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
hold off
if savePlots
    if useDefault
        saveas(h,[StatsPath tagfolder filesep 'modesProfileAllSelected' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
    else
        saveas(h,[StatsPath tagfolder filesep 'modesProfileAll' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXPORT MODES CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
add2selection=1
if cluster2movit
    nameSel_modes=[datapath dataset tagfolder '_t' filesep 'modesClusters' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) xtag '.csv']
    fid=fopen(nameSel_modes,'w');
    fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s\n','id_center','selection','timestep','x','y','z','id_mother','validation');
    for cl=1:length(classes)
        Class=ClusterData(C==cl,:,:);
        'Class size'
        cl
        sclass=size(Class,1)
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

% /* --------------------------------------------------------------------------------------
%  * File:    runModes.m
%  * Date:    01/12/2017
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2013-2017, David Pastor Escuredo
%  with Biomedical Image Technology, UPM (BIT-UPM)
%  with BioEmergences, CNRS
%  with LifeD lab
%  All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering profiles together %%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Remove invalid outliars in a material perspective that will affect the clustering
            %Filtering to smooth if necessary
            toremove=[];
            c=1;
            numOfMaterialProfiles=size(clus)
            %general statistical cropping
            mx=prctile(clus(:),98)
            mn=prctile(clus(:),2)
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
            
            save([desc 'profile-' tagDesc '-' selCtag], 'clus');
            
            
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        else
            load([desc 'Basis' tagDesc seltagBasis tag_modes '-modes' num2str(mxCluster)])
        end
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
            rangeD=[0 prctile(clus_old(:),90)];
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
                colorrgb=modes_colors{DescriptorIndex}(icc,:)/255
                plot(Basis(icc,:)','LineWidth',widthMode,'Color',colorrgb);
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
        for ib=1:size(Basis,1)
            mm_lg{ib}=[nameDescLeg '_{m' num2str(ib) '}'];
        end
        if mxCluster>3
        legend(mm_lg,'Location','NorthOutside','Orientation','horizontal','FontSize',fsize-14)
        else
        legend(mm_lg,'Location','NorthOutside','Orientation','horizontal','FontSize',fsize)
            
        end
        %Save
        screen_size = get(0, 'ScreenSize');
        set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        hold off
        if savePlots
            if useDefault
                saveas(h,[StatsPath dataset tagfolder filesep 'modesProfileSelected' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
            else
                saveas(h,[StatsPath dataset tagfolder filesep 'modesProfile' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
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
end

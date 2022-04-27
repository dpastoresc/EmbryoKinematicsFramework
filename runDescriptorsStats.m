% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stastistics package run script

clear;
StatsPath='../MechanicsStats/'
v_tag=''

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User's Configuration %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% List of datasets to process. datasetListBioEmergences
dNs=[1]%4ataset index to load info
tinilags=[51 62 39 8 39]%time refernece for cumulative deformation analysis

%%%%%%%%%%%%%%%% Descriptor Set and Descriptor Index %%%%%%%%%%%%%%%%%%%%%
%1 Topology
%2 Strain rates
%3 Displacements Field descriptors (Unstable)
%4 Left Lagrangian descriptors
%5 Right Lagrangian descriptors
DescriptorsSets=5
DescriptorsIndexes=[5]

% For 4 and 5:
%{'FTLE','FTLE isocoric','Speed Lagrangian','Topology','J','MIC1','MIC2','Rotation','Shear Angle Abs','Shear Angle Intermediate','e1','e2','e3'}
if DescriptorsSets>3
DescriptorsIndexes=[5 6 7 8]
elseif DescriptorsSets==1
% For 1:
%{'Speed','Speed Averaged','Topology','P','Q','D','Expansion planar', 'Compression planar', 'Rotation rate'}
DescriptorsIndexes=[4 5 6 ]
elseif DescriptorsSets==2
% For 2:
%{'Qs','Qd','Residue Qd','Max Shear Angle','e1','e2','e3','d1','d2','d3'}
DescriptorsIndexes=[1 2]
end

%%%%%%%%%%%%%%% Descriptor Kinematics. loadStatsMetadata
%%%%%%%%%%%%%%%% Running options
%Do Domain profiling
doProfile=1
old=0
xtag='-1-0'
drawProfile=1
config=222%Configuration of selections colors and so on. Check selection options
simple=2%0- all detail 1- simple 2- nothing
min_samples=30

%Mono descriptors clustering of Domain
doCluster=0
drawClusterProfile=0

%%%%%%%%%%%%%%%% Tissue/Cell selection
%Tissue selection
withSelection=1
seltag='-epi'
extraTag=''
seltag=[seltag extraTag]
selections=[1 2 3 4]
selections=[1 2 6 4]
%selections=[1 3]

%%%%%%%%%%%%%%%% New limits in time
%-1 is default interval of the datasets
tfinal_m=-1
tini_m=-1
hini=6.8
hfin=14
hpf_limits=1
crop_time_axis=1
phys_time=1
plotZero=1
plotHPF=1%Flag
t_max_cluster=-1%178
printTrefs=1
startInTref=1
setYaxis=0
writeZero=1

%%%%%%%%%%%%%%%% Stats Options PROFILE
make_hist=0%flag to get the histogram of each selection descriptor
spatial_mask=1%1-all ids in each time (USE THIS) %0- we rely on material ids without considerint time
tref_propagate=-1%we generate a material mask using the ids of a specific timestep
%Statistical cropping for each selection
removeOutliars=0;%1-remove 0-saturate. perenctiles are configure in statsMetadata for each descriptor
% Temporal gradient of spatial information
derive=0%flag
spatialRad=15%tolerance

%%%%%%%%%%%%%%% Stas Options DRAW PROFILE
fsize=40
gridOn=0
gridStyle='--'
hidePlot=0
savePlot=1
downsampleT=0
overlayed=0

%%%%%%%%%%%%%%%% Stats Options CLUSTER Mod
%cluster over specific selections
only_selection=withSelection%considers only the selection
clusterSelec=[3]
%Hierarchical clustering
distanceT='Euclidean';
method='average';
mxCluster=6%length(clusterSelec)*4;
saveClusterData=0
%filtering
removeProfileWithInvalid=1
smoothN=11%window for 1D filtering of profiles
%plots
seePlots=1
savePlots=1
%export to movit selections back
cluster2movit=1
add2selection=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doBasis=0
doProjection=0
useDefault=0
if startInTref
    phys_time=0
    plotHPF=0
end
for idN=1:length(dNs)   
    'Dataset' 
    dN=dNs(idN)
    tinilag=tinilags(dN)
    dataset
    loadParametersKinematics;
    for DescriptorsSet=DescriptorsSets
        for DescriptorIndex=DescriptorsIndexes            
            %nameSel=[datapath dataset '_t_sim/' dd '_t' seltag '_all_T']
            nameSel=[datapath dataset '_t/' dd '_t' seltag '_all_T']     
            debug_this=0
            invalidValue=-1
            loadStatsMetadata
            if withSelection
                selec=dlmread([nameSel '.csv'], ';', 1, 0);
            else
                seltag='-all'
            end
            if (doCluster && cluster2movit) || drawClusterProfile
                nameSel_clus=[nameSel '-' selCtag tagsCmaps '-tf' num2str(t_max_cluster) '-s' num2str(mxCluster) '.csv']
            end              
           
            if doProfile
                if old
                    runDescriptorsProfile2
                else
                    metric=[]          
                    for si=1:length(selections)
                       ss=selections(si)
                       runDescriptorsProfile
                    end
                    save(stats_file, 'metric');
                end
            end
          
            if drawProfile
                drawDescriptorsProfile
            end
            if doCluster
                runDescriptorsCluster
            end
            if drawClusterProfile
                config=99
                %loadParametersKinematics;
                selec=dlmread([nameSel_clus], ';', 1, 0);
                drawDescriptorsProfile
            end       
        end
    end
end




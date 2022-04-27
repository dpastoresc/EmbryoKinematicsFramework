% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stastistics package run script
clear;
StatsPath='../MechanicsStats/'
v_tag=''

modes_colors={[97 97 97; 174 119 40; 150 150 150;],[97 97 97; 174 119 40; 150 150 150;],{},{},[97 97 97; 174 119 40; 150 150 150;],[0 0 143; 218 179 255; 191 0 191],[20 20 83; 238 209 175; 211 0 121],[0 153 153; 0 204 0; 0 255 100],[97 97 97; 174 119 40; 150 150 150;],[97 97 97; 174 119 40; 150 150 150;],[97 97 97; 174 119 40; 150 150 150;],[97 97 97; 174 119 40; 150 150 150;],[97 97 97; 174 119 40; 150 150 150;],[97 97 97; 174 119 40; 150 150 150;]}
%%%%%%%%%%%%%%%%%%%%%%%%%% User's Configuration %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% List of datasets to process. datasetListBioEmergences
dNs=[1]
tinilags=[61 62 39 8 39 95]%for each dataset
tinilags(13)=50
tinilags(18)=69
tinilags(14)=105
%t_max_cluster=250
%t_max_cluster=213
%t_max_cluster=194
redo=1
%071226a
%61 and 118 155 -> 178
%081018a
%95 and 154 190-> 213
%091021aF
%105 154 -> 204

%%%%%%%%%%%%%%%% Descriptor Set and Descriptor Index %%%%%%%%%%%%%%%%%%%%%
%1 Topology
%2 Strain rates
%3 Displacements Field descriptors (Unstable)
%4 Left Lagrangian descriptors
%5 Right Lagrangian descriptors
DescriptorsSets=[1]
%For 4-5
DescriptorsIndexes=[5 6 7 8]%check libStats/loadStatsMetadata
%For 1
DescriptorsIndexes=[4]%check libStats/loadStatsMetadata
%For 2
%DescriptorsIndexes=[1 2]%check libStats/loadStatsMetadata

%%%%%%%%%%%%%%%% Tissue selection from Movit. Original domain to get modes
%Tissue selection to get the modes
withSelection=1%keep 1
seltagBasis='-all'%all71'%Selection of the material domain basis
clusterSelecBasis=[4]%Selec numbers
xtag='-1-0'
old=0
reclimit=480

%%%%%%%%%%%%%%%% Configuration Modes by Unsupervised Clustering
useDefault=1;%Number of modes in loadStatsMetadata
mxCluster=4;%Number of modes to find
clus_method=0%0: K-means uses distanceK 1: Hierarchical clustering uses distanceT
distanceType=2%0:correlation 1: euclidean 2: cosine 3: mahalanobis
smoothN=5%window for 1D filtering of profiles before clustering (timesteps)
cluster2movit=0%export clustering after monomodal or multimodal segmentation into a new label file for MovIt
exportBasisMap=0%Distances map when generating the modes to Movit
%Kmeans options
startCluster=0%if 1 it uses 10% of samples for preliminary clustering
%Hierarchical options
methodHier='ward'
tagDis=''
if distanceType==1
    distanceK='sqEuclidean'
    distanceT='Euclidean';
    tagDis='euc'
elseif distanceType==0
    sameRangeDis=1
    distanceK='correlation'
    distanceT='correlation'
    tagDis='cor'
elseif distanceType==2
    distanceK='cosine'
    distanceT='cosine'
    tagDis='cos'
elseif distanceType==3
    distanceK='mahalanobis'
    distanceT='mahalanobis'
    tagDis='mah'
end

%%%%%%%%%%%%%%%%% Plotting options
extraPlots=0
seePlots=1
savePlots=1
textSimple=1
fsize=36

%%%%%%%%%%%%%%%%% time
tfinal_m=-1
tini_m=-1
hini=6
hfin=13
hpf_limits=1
crop_time_axis=1
phys_time=1
plotZero=0
plotHPF=1%Flag

step_time=20
gridOn=1
widthMode=7
use_defaultRange=1
saveClusterData=0%save
removeProfileWithInvalid=1
%%%%%%%%%%%%%%%%%%%%% User's Configuration End %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% RUNNING %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selections=[1 2 3 1]
only_selection=withSelection%considers only the selection
derive=0
doCluster=0
doProfile=0
doBasis=1
doProjection=0
drawClusterProfile=0

for idN=1:length(dNs)
    'Dataset'
    dN=dNs(idN)
    tinilag=tinilags(dN)
    loadParametersKinematics;
    for DescriptorsSet=DescriptorsSets
        %for DescriptorIndex=DescriptorsIndexes
        %nameSel=[datapath dataset '_t_sim/' dd '_t' seltag '_all_T']
        
        debug_this=0
        invalidValue=-1
        loadStatsMetadata
        
        nameSel=[datapath dataset '_t/' dd '_t' seltagBasis '_all_T']
        if withSelection
            selec=dlmread([nameSel '.csv'], ';', 1, 0);
        else
            seltag='-all'
        end
        selCtag=selection2tissue(clusterSelecBasis,'');
        runModes
        
    end
end



    
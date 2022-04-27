% /* --------------------------------------------------------------------------------------
%  * File:    runMechanomeSemengtation.m
%  * Date:    01/12/2017
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2013-2018, David Pastor Escuredo
%  with Biomedical Image Technology, UPM (BIT-UPM)
%  with BioEmergences, CNRS
%  with LifeD lab
%  All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Stastistics package run script
clear;
StatsPath='../MechanicsStats/'
v_tag=''

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User's Configuration %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% List of datasets to process. datasetListBioEmergences
dNs=[1]
tinilags=[61 60 39 8 39 190]%for each dataset
tinilags(13)=50
tinilags(18)=69
% t_max_cluster=178
% t_max_cluster=213
% t_max_cluster=194
reclimit=1100
limitP=0

%%%%%%%%%%%%%%%% Descriptor Set and Descriptor Index %%%%%%%%%%%%%%%%%%%%
%1 Topology
%2 Strain rates
%3 Displacements Field descriptors (Unstable)
%4 Left Lagrangian descriptors
%5 Right Lagrangian descriptors
DescriptorsSets=[5]
%For 4-5
DescriptorsIndexes=[5 6 8]%6 7 8]%check libStats/loadStatsMetadata
%For 1
%DescriptorsIndexes=[4 5 6 9 ]%check libStats/loadStatsMetadata
%For 2
%DescriptorsIndexes=[1 2]%check libStats/loadStatsMetadata

%%%%%%%%%%%%%%% Segmentation parameters
%Distance of each profile to the mode
distanceClass='euclidean'
%Do Segmetnation or just projection
segmentField=1%Segments the projected selection using the modes and a hierarchical clustering
%Number of segmentation
segmentationClasses=4%number of clusters for segmentation of the field projected
%Export segmentation to Movits
cluster2movit=1
%Project by amount of tissue -DEPRECATED
granulometry=1
%Options to work with the distances map
cropDisRange=1
sameRangeDis=1
binarizeProjection=1
xtag='-1-0'
old=0
transp=0

%%%%%%%%%%%%%%%% Extra
method='ward';%Hierarchical clustering for segmentField=1
saveClusterData=0%save
removeProfileWithInvalid=1

%%%%%%%%%%%%%%%% Tissue selection from Movit
%Tissue selection to get the modes
withSelection=1%keep 1
seltagBasis='-all'%all71'%Selection of the material domain basis
clusterSelecBasis=[4]%Selec numbers

%Tissue selection to project
seltag='-all'%Selection of the material domain to project
clusterSelec=[4]%Selec
extraTag=''
seltag=[seltag extraTag]

%%%%%%%%%%%%%%%%% Plotting options
extraPlots=0
seePlots=0
savePlots=1
textSimple=1
fsize=36

%%%%%%%%%%%%%%%%% Time
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

%%%%%%%%%%%%%%%% Configuration of the Modes Used.
useDefault=1;%Number of modes in loadStatsMetadata
mxCluster=-1;%Number of modes to find
clus_method=0%0: K-means uses distanceK 1: Hierarchical clustering uses distanceT
distanceType=2% 0:correlation 1: euclidean
tagDis=''
if distanceType==1
    distanceK='sqEuclidean'
    distanceT='Euclidean';
    tagDis='euc'
elseif distanceType==0
    distanceK='correlation'
    distanceT='correlation'
    tagDis='cor'
elseif distanceType==2
    distanceK='cosine'
    distanceT='cosine'
    tagDis='cos'
end
smoothN=1%window for 1D filtering of profiles before clustering (timesteps)

%%%%%%%%%%%%%%%%%%%%% User's Configuration End %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% RUNNING %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modes_colors={{},{},{},{},[97 97 97; 174 119 40; 150 150 150],[0 0 143; 218 179 255; 191 0 191],[20 20 83; 238 209 175; 211 0 121],[0 153 153; 0 204 0; 0 255 100]}
selections=[1 2 3 4]
only_selection=withSelection%considers only the selection
derive=0
doCluster=0
doBasis=0
doProjection=1
doProfile=0
drawClusterProfile=0

%
for idN=1:length(dNs)
    'Dataset'
    if old
        xtag=''
    end
    dN=dNs(idN)
    tinilag=tinilags(dN)
    dataset
    loadParametersKinematics;
    for DescriptorsSet=DescriptorsSets
        %for DescriptorIndex=DescriptorsIndexes
        %nameSel=[datapath dataset '_t_sim/' dd '_t' seltag '_all_T']
      
        debug_this=0
        invalidValue=-1
        loadStatsMetadata
        DescriptorsSets
        dNs
        length(dNs)
          
        nameSel=[datapath dataset '_t/' dd '_t' seltag '_all_T']
        if withSelection
            selec=dlmread([nameSel '.csv'], ';', 1, 0);
        else
            seltag='-all'
        end
        selCtag=selection2tissue(clusterSelec,extraTag);
        runMechanomeDomain
    end
end




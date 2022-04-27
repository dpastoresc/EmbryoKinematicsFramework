% /* --------------------------------------------------------------------------------------
%  * File:    runTrackingTensors.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

clear
%%%%%%%%%%%%%%%%%%%%% SELECT DATASET %%%%%%%%%%%%%%%%%%%%
dN=1%datasetListBioEmergences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tinilag=-1%not used here
xtag=''
loadParametersKinematics;


%Tracking format
%0 EMB FILE FORMAT PROCESSING
if inputTrackingType==0
    
    datasetBioemergences= dataset%[dataset '_tn' trackID]
    embfilename = [datasetBioemergences 't' trackID '.emb']
    
    build_mov_F_E_R = 1;
    Rscale=40%not used but needs to be init
    radialCorrection=0
    
    generateArrays=0;
    if redoDataStructure
        generateArrays=1;
    else
        if ~exist([datapath dataset '_t' filesep datasetBioemergences '_XYZ.mat' ])
            'Need to create structures'
            generateArrays=1;
        end
    end
    
    %Very specific EMB format flags (for version compatibility)
    indexZero=1;%some data started with timestep=1, but new data starts with 0. carefull this requires some
    numCol=9;
    zoff=0;
    %Adapt to matlab indexes
    tiniMat=tini+1
    tfinalMat=tfinal+1

    if generateArrays
        'Create tracking data structures from tracking EMB file'
        emb2dataArrays;
    end
    
    %computation of continuum and discrete tensors
    %this works on the data arrays in datapath  
    launch_EMBlib
end


'Tensors calculated'
runEulerianDescriptors
'Eulerian descriptors calculated'

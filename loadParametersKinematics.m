% /* --------------------------------------------------------------------------------------
%  * File:    loadParametersKinematics.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

%Welcome!
%Set up the configuration for getting kinematics descriptors of your tracking data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER  CONFIGURATION %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SET!! Tracking DATA
%For now we support BioEmergences EMB tracking format
%0->emb BIOEMERGENCES
%>0->other formats to implement
inputTrackingType=0;   
if inputTrackingType==0
    % Emb files
    %       cellid;cellnum;x;y;z;t;motherref;validation
    if dN<0
    dN=1%THIS IS THE DATASET POINTED OUT IN datasetList.m
    end
    export2movit=1;
    loadDatasetBioEmergences;
end

%Configure loadParametersTrackingTensors for low level configuration of the
%tracking to continuum processing
loadParametersTrackingTensors

%%%%%%%%%%%%%%%%%%%%%%%%%% SET!!Langrangian parameters to use if langrange=1
%Integration window in minutes
%if 0, the window increments along time
%if >0, the window is fix
integrationT=0;
if tinilag<0
    tinilag=51;%the starting point to accumulate the data
end
%Time integration process option
update_pos=1%Update the prediction of the sample to observe the closest future cell to link
weight_similar_vel=0;%considers the velocity vector additionally to the distance to link trajs.
save_tracking=1;%Save the dynamic lagrangian trackings from tinilag to tfinal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END USER CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% Time configuration %%%%%%%%%%%%%%%%%%%%
%You can specify tini and tfinal in loadDataset, otherwise will take the
%limits of the tracking data
if tini==-1
    tini=0;
end
if tfinal==-1
    tfinal=steps-1;
end
step=1;
vsteps=tfinal-tini+1;
T=uint16(integrationT*60*1000/dt_ms)%Integration window in time steps
T=double(T);%otherwise we break it all

%%%%%%%%%%%%%%%%%% Configuration of paths. Do not touch
Xave=Xscale
Tave=Tscale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All naming here
'Configuring Done: Creating Paths...'
d=dataset;
%depends on how the tensors are calculated.s
if onlySimilar
    folder =[ datapath dataset '_tn' trackID '_t/']
else
    folder =[ datapath dataset '_tn' trackID '_t_test/' ]
end
if dataLocal
    folder =[ datapath dataset '_t/']
end

folderkin=[folder '/kinematics_' dataset filesep]
desc=[folder '/kinematics_' dataset filesep]
desc1=folder
%desc=folder%check this!!!
mkdir(desc);
dataset=[d '_tn' trackID];
if dataLocal
    dataset=d
end
param=['T' num2str(Tave) '_X' num2str(Xave) '_G' num2str(X)];
if ~spatialAve
    tensorid=['F_' param];
    strainid=['E_' param];
else
    tensorid=['Fave_' param];
    strainid=['Eave_' param];
end
neighid=['Neigh_' param];
veloid=['Mov_' param];
velomeanid=['MeanMov_' param];
veloRid=['RMov_' param];
velomeanslid=['MeanMovSl_' param];
veloslid=['MovSl_' param];

tagvel='';
if spatialAve
    tagvel='_ave';
end
rawTopoDescriptors=[desc dataset tensorid '_topo.mat'];
rawStrainDescriptors=[desc dataset strainid '_strains.mat'];
rawVelDescriptors=[desc dataset veloid '_vel.mat'];
tags=[''];
if T==0
    tags=['_' num2str(tinilag)];
end
if weight_similar_vel
    rawLagDescriptors=[desc dataset tagvel param '_lag_T' num2str(T) tags '_sim' xtag '.mat'];
    rawLagDescriptorsLeft=[desc dataset tagvel param '_lagleft_T' num2str(T) tags '_sim' xtag '.mat'];
else
    rawLagDescriptors=[desc dataset tagvel param '_lag_T' num2str(T) tags xtag '.mat'];
    rawLagDescriptorsLeft=[desc dataset tagvel param '_lagleft_T' num2str(T) tags xtag '.mat'];
end


% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%%%%%%%%%%%%%%%%%%%%%% SELECT YOUR DATASET %%%%%%%%%%%%%%%%%
%USER!!!!
%Set your Dataset id in datasetList.m
%dN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dN=1
%tinilag=-1
loadLibs;
loadDatapath;
datasetListBioEmergences;

%%% Presets
dd=dataset
sizePixel=[512 512 512]
spacingPixel=[1 1 1]
steps=100;
vsteps=steps;
t0=0;%secs
dt_ms=1;%1 milisec
species='undefined'
bioemergences=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% BIOEMERGECES DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%DatasetInfo%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sizePixel      pixels in x y z
%spacingPixel   spacing in microns for x y z
%steps          number of timesteps in the time-lapse
%t0             start of imaging in "
%dt             timestep in milisecs betweent steps
%species        animal of the dataset 
if bioemergences && ~dataLocal
    [sizePixel spacingPixel steps t0 dt_ms species ]=readDataParamBioemergences(dataset);
    if tini==-1
        tini=0;
    end
    if tfinal==-1
        tfinal=steps-1;
    end
    if steps<(tfinal+1)
        tfinal=steps-1;
    end
    vsteps=tfinal+1;
end

t0_ms=t0*1000;
dt_hours=dt_ms/(1000*3600);%from milisec to hours
t0_hours=t0_ms/(1000*3600);%from miliec to hours
dt_min=dt_ms/(1000*60)%from milisec to min
t0_min=t0_ms/(1000*60);%from miliec to min
dt_sec=dt_ms/(1000);%from milisec to min
t0_sec=t0_ms/(1000);%from miliec to min



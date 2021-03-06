% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved


%%%% Stats metadata
stat={'Mean','Median','std','Max',',Min','Percentile-99','Percentile-1' ...
    'Percentile-75','Percentile-25','entries','entries-dif','samplesTaken' ...
    'samplesInside', 'ss', 't'} 
   
tagfolder=''
%Topology Descriptors
tagg=''
if DescriptorsSet==1
    tags={'Speed','Speed Averaged','Topology','P','Q','D','Expansion planar', 'Compression planar', 'Rotation rate'}
    tagsShort={'Speed','SpeedAve','Topo','P','Q','\delta\alpha','ExpPlanar','CompPlanar', 'RotRate'}
    tagsCmaps=tags;
    index=[7 8 12 13 14 15 16 20 24]
    if doCluster || doProfile || doBasis || doProjection
        load(rawTopoDescriptors); %'DescriptorsT');
        Desc=DescriptorsT;
        clear DescriptorsTS;
    end
    desc_limits=[-0.0005 0.0005; -0.0005 0.0005; -1 -1; -0.01 0.01; -0.00012 0.00012; -0.000000000001 0.000000000001; -1 -1; -1 -1; 0 1]
    pmax=[95 95 95 85 85 80 95 95 95]
    pmin=[5  5  5  15 15 20 5  5  5]
    zv=[0 0 0 0 0 0 0 0]
    isSymmetric=[0 0 0 1 1 1 0 0 0]
    descriptorsIndexes=[4 5 6 9 ]%check libStats/loadStatsMetadata
    
    climits=[-1 -1; -1 -1; -1 -1; -1 -1; -0.2 0.2; 0 1; 0 1; 0 30; -1 -1; -1 -1; -1 -1; -1 -1; -1 -1]
    modes_mech=[0 0 0 0 3 3 3 3 3 3 3 3 3]

elseif DescriptorsSet==2
    tags={'Qs','Qd','Residue Qd','Max Shear Angle','e1','e2','e3','d1','d2','d3'}
    tagsShort={'Qs','Qd','Res','MaxShear','e1','e2','e3','d1','d2','d3'}
    tagsCmaps=tags;
    index=[7 20 22 21 8 12 16 23 27 31]
    if doCluster || doProfile
        load(rawStrainDescriptors)%, 'DescriptorsS')
        Desc=DescriptorsS;
        clear DescriptorsS;
    end
    desc_limits=[0 0.0001; 0 0.0001; -1 -1; -0.01 0.01; -0.00015 0.00015; -0.000000000005 0.000000000005; -1 -1; -1 -1; 0 1]
    pmax=[95 85 95 95 85 80 95 95 95]
    pmin=[5  15  5  5 15 25 5  5  5]
    zv=[0 0 0 0 0 0 0 0]
    isSymmetric=[0 0 0 0 1 1 1 1 1 1]
    DescriptorsIndexes=[1 2]%check libStats/loadStatsMetadata

elseif DescriptorsSet==3
    tags={'Qs','Qd','Residue Qd','Max Shear Angle','e1','e2','e3','d1','d2','d3'}
    tagsShort={'Qs','Qd','Res','MaxShear','e1','e2','e3','d1','d2','d3'}
    index=[7 20 22 21 8 12 16 23 27 31]
    if doCluster || doProfile
        load(rawStrainDescriptors)%, 'DescriptorsS')
        Desc=DescriptorsV;
        clear DescriptorsV;
    end
    desc_limits=[-0.0005 0.0005; -0.0005 0.0005; -1 -1; -0.015 0.015; -0.0002 0.0002; -0.000000000001 0.000000000001; -1 -1; -1 -1; 0 1]
    pmax=[95 95 95 95 85 80 95 95 95]
    pmin=[5 5 5 5 15 20 5 5 5]
    zv=[0 0 0 0 0 0 0 0]
    isSymmetric=[0 0 0 0 1 0 0 0]
    DescriptorsIndexes=[1 2]%check libStats/loadStatsMetadata
    
elseif DescriptorsSet==4
    tagsCmaps={'FTLE','FTLEisochoric','Speed Lagrangian','Topology','J','MIC1','MIC2','Rotation','MaxShear','IntermediateShear','e1','e2','e3'}
    tags={'FTLE','FTLEisochoric','Speed Lagrangian','Topology','\Delta V','\Delta sigma1','\Delta sigma2','\theta','MaxShear','IntermediateShear','e1','e2','e3'}
    if useDefault
        tagsShort={'FTLE','FTLEiso','SpeedLag','Topo','\Delta V','\Delta sigma1','\Delta sigma2','\theta','Shear1','Shear2','e1','e2','e3'}
    else
        tagsShort={'FTLE','FTLEiso','SpeedLag','Topo','\Delta V','\Delta sigma1','\Delta sigma2','\theta','S1','S2','e1','e2','e3'}
    end
    %index=[7 8 9 10 11 12 13 29 17 21 25]
    index=[7 8 9 10 11 12 13 14 30 31 18 22 26]
    if doCluster || doProfile
        load(rawLagDescriptors)%, 'DescriptorsV'
        Desc=DescriptorsL;
        clear DescriptorsL;
    end
    if T>0
        tagg=[tagg '-lag-' num2str(T)]
    else
        tagg=[tagg '-lag-t' num2str(tinilags(dN))]
    end
     isSymmetric=[0 0 0 0 1 0 0 0 0 0 1 1 1]
    zv=[0 0 0 0 1 0 0 0 0 0 1 1 1 ]
    pmax=[95 95 95 100 100 90 90 90 90 90 90  90 10 10]
    pmin=[5  5  0   0  0 0   0  0  0  0  0  10 10 10]
    %colormap is consistent with desc_limits
    climits=[-1 -1; -1 -1; -1 -1; -1 -1; -0.2 0.2; 0 1; 0 1; 0 30; -1 -1; -1 -1; -1 -1; -1 -1; -1 -1]
    modes_mech=[0 0 0 0 3 3 3 3 3 3 3 3 3]

elseif DescriptorsSet==5
    %old tags={'FTLE','FTLE isocoric','FTLE backwards','FTLE isocoric backwards','Topology','J','Rotation','Distortion','e1','e2','e3'}
    %tags={'FTLE','FTLE isocoric','Speed Lagrangian','Topology','J','MIC1','MIC2','\theta','MaxShear','IntermediateShear','e1','e2','e3'}
    tagsCmaps={'FTLE','FTLEisochoric','Speed Lagrangian','Topology','J','MIC1','MIC2','Rotation','MaxShear','IntermediateShear','e1','e2','e3'}
    tags={'FTLE','FTLEisochoric','Speed Lagrangian','Topology','\Delta V','\Delta \gamma _1','\Delta \gamma _2','\Delta \alpha','MaxShear','IntermediateShear','e_1','e_2','e_3'}
    if useDefault
        tagsShort={'FTLE','FTLEiso','SpeedLag','Topo','\Delta V','\Delta \gamma _1','\Delta \gamma _2','\Delta \alpha','Shear1','Shear2','e_1','e1_2','e_3'}
    else
        tagsShort={'FTLE','FTLEiso','SpeedLag','Topo','\Delta V','\Delta \gamma _1','\Delta \gamma _2','\Delta \alpha','S1','S2','e_1','e_2','e_3'}
    end
    index=[7 8 9 10 11 12 13 14 30 31 18 22 26]
    if doCluster || doProfile || doBasis || doProjection
        load(rawLagDescriptorsLeft)
        Desc=DescriptorsLLeft;
        clear DescriptorsLLeft;
    end
    if T>0
        tagg=[tagg '-lagL-' num2str(T)]
    else
        tagg=[tagg '-lagL-t' num2str(tinilag)]
    end
    desc_limits=[-1 -1; -1 -1; -1 -1; -1 -1; 0.7 1.2; 0 1.3; 0 1.3; 0 35]
    desc_limits=[-0.05 0.05; -0.05 0.05; -1 -1; -1 -1; 0.7 1.2; 0 1.6; 0 1.6; 0 50; 0 0.6; 0 0.6; -1 -1; -1 -1; -1 -1]
   
    isSymmetric=[0 0 0 0 1 0 0 0 0 0 1 1 1]
    zv=[0 0 0 0 1 0 0 0 0 0 1 1 1 ]
    pmax=[95 95 95 100 80 90 90 90 90 90 90  90 10 10]
    pmin=[5  5  0   0  20 0   0  0  0  0  0  10 10 10]
    %colormap is consistent with desc_limits
    climits=[-1 -1; -1 -1; -1 -1; -1 -1; -0.2 0.2; 0 1; 0 1; 0 30; -1 -1; -1 -1; -1 -1; -1 -1; -1 -1]
    modes_mech=[0 0 0 0 3 3 3 3 3 3 3 3 3]
    
    if DescriptorsIndexes(1)==1
      DescriptorsIndexes=[5 6 7 8 9 10 1 2]%check libStats/loadStatsMetadata
    end
    %DescriptorsIndexes=[5]
end

tagmain='profile'
if derive
    tagmain='gradientSpatialInfo'
else
    spatialRad=0
end

if tfinal_m>-1
    t_limit=tfinal_m
else
    t_limit=tfinal
end
if tini_m>-1
    t_start=tini_m
else
    t_start=tini
end
%Lag->configure time limits
if DescriptorsSet>3
    %lag
    if T==0
        t_start=tinilag
    else
        t_T=tfinal-T+1
        if t_limit>t_T
            t_limit=t_T
        end
    end
end


%%%%%%%%%%%%%%%%% Names %%%%%%%%%%%%%%%%
mkdir([StatsPath dataset]);
if ~doBasis && ~doProjection
    tagDescOr= tagsCmaps{DescriptorIndex}
    tagDesc=['-' tagsCmaps{DescriptorIndex} tagg]
    descriptor=index(DescriptorIndex)
    stats_options=[tagmain tagDesc seltag v_tag '-' num2str(pmax(DescriptorIndex)) '-' num2str(pmin(DescriptorIndex))]
    stats_file=[StatsPath dataset tagfolder filesep stats_options]
    profile_file=[StatsPath stats_options '-' dataset]
end

if drawClusterProfile
    profile_file=[profile_file '-cluster']
end

%%%% HPF bars
hpf6=6*60
hpf7=7*60
hpf8=8*60
hpf9=9*60
hpf10=10*60
hpf11=11*60
hpf12=12*60
hpf13=13*60
hpf14=14*60
hl=1.5
hc=[0.7 0.4 0.95]


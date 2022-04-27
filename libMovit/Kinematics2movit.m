% /* --------------------------------------------------------------------------------------
%  * File:    Kinematics2movit.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

%Sets up the data to create Movit Descriptors files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Dataset parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loadParametersKinematics;
loadParametersKinMovit;%import parameters to run the translation to movit

%Insert colormaps
color_folder='colormaps/'
cdir=dir(color_folder);

%%%%%%%%%%%%%%%% SubPaths %%%%%%%%%%%%%%%%%%%%%%%%%%%
if vheader==2
    desc=[desc]%'v2/'];
    mkdir(desc);
end
if topology
    if spatialAve
        descF=[desc 'Fave_' param];
    else
        descF=[desc 'F_' param];
    end
    if globalNorm
        descF=[descF '_norm'];
    end
    mkdir(descF);
    fileHeader=[descF '/header.kdh']
    tlapse=[descF '/timelapse/'];
    mkdir(tlapse);
    if vheader==2
        fileHmax=[descF '/maxtemp.kdh']
        fileHmin=[descF '/mintemp.kdh']
    end
    for ii=3:length(cdir)
        if length(strfind(cdir(ii).name,'.png'))>0
            copyfile([color_folder cdir(ii).name], descF);
        end
    end
    TopologyDescriptors2movit
end
if strains
    if spatialAve
        descE=[desc 'Eave_' param];
    else
        descE=[desc 'E_' param];
    end
    if globalNorm
        descE=[descE '_norm'];
    end
    mkdir(descE);
    fileHeader=[descE '/header.kdh']
    tlapse=[descE '/timelapse/'];
    mkdir(tlapse);
    if vheader==2
        fileHmax=[descE '/maxtemp.kdh']
        fileHmin=[descE '/mintemp.kdh']
    end
    
    tlapse=[descE '/timelapse/'];
    mkdir(tlapse);
    for ii=3:length(cdir)
        if length(strfind(cdir(ii).name,'.png'))>0
            copyfile([color_folder cdir(ii).name], descE);
        end
    end
    StrainsDescriptors2movit
end
if velocity
    descV=[desc 'V_' param];
    if globalNorm
        descV=[descV '_norm'];
    end
    mkdir(descV);
    fileHeader=[descV '/header.kdh']
    tlapse=[descV '/timelapse/'];
    mkdir(tlapse);
    if vheader==2
        fileHmax=[descV '/maxtemp.kdh']
        fileHmin=[descV '/mintemp.kdh']
    end
    tlapse=[descV '/timelapse/'];
    mkdir(tlapse);
    for ii=3:length(cdir)
        if length(strfind(cdir(ii).name,'.png'))>0
            copyfile([color_folder cdir(ii).name], descV);
        end
    end
    MigDescriptors2movit
end


%DEPRECATED
if lagrange
    left=0
    descL=[desc 'L_' param '_w' num2str(T)];
    mkdir(descL);
    fileHeader=[descL '/header.kdh']
    tlapse=[descL '/timelapse/'];
    mkdir(tlapse);
    if vheader==2
        fileHmax=[descL '/maxtemp.kdh']
        fileHmin=[descL '/mintemp.kdh']
    end
    tlapse=[descL '/timelapse/'];
    mkdir(tlapse);
    lDescriptors=rawLagDescriptors
    LagrangeDescriptors2movit
    left=1
    descL=[desc 'C_' param '_w' num2str(T)];
    mkdir(descL);
    fileHeader=[descL '/header.kdh']
    tlapse=[descL '/timelapse/'];
    mkdir(tlapse);
    if vheader==2
        fileHmax=[descL '/maxtemp.kdh']
        fileHmin=[descL '/mintemp.kdh']
    end
    tlapse=[descL '/timelapse/'];
    mkdir(tlapse);
    for ii=3:length(cdir)
        if length(strfind(cdir(ii).name,'.png'))>0
            copyfile([color_folder cdir(ii).name], descL);
        end
    end
    lDescriptors=rawLagDescriptorsLeft
    LagrangeDescriptors2movit
end



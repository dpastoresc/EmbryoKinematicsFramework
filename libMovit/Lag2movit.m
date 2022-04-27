% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%Sets up the data to create Movit Descriptors files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Dataset parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loadParametersKinematics;
loadParametersKinMovit;%import parameters to run the translation to movit

%%%%%%%%%%%%%%%% SubPaths %%%%%%%%%%%%%%%%%%%%%%%%%%%
if vheader==2
    desc=[desc]%'v2/'];
    mkdir(desc);
end

timelapse='/timelapse/'

if FixedWindow    
    left=0
    descL=[desc 'L_' param '_w' num2str(T) '-' num2str(real_pos) '-' num2str(reuse_links)];
    mkdir(descL);
    fileHeader=[descL '/header.kdh']
    tlapse=[descL timelapse];
    mkdir(tlapse);
    if vheader==2
        fileHmax=[descL '/maxtemp.kdh']
        fileHmin=[descL '/mintemp.kdh']
    end
    tlapse=[descL timelapse];
    mkdir(tlapse);
    lDescriptors=rawLagDescriptors
    LagrangeDescriptors2movit  
    left=1
    descL=[desc 'C_' param '_w' num2str(T) '-' num2str(real_pos) '-' num2str(reuse_links)];
    
    mkdir(descL);
    fileHeader=[descL '/header.kdh']
    tlapse=[descL timelapse];
    mkdir(tlapse);
    if vheader==2
        fileHmax=[descL '/maxtemp.kdh']
        fileHmin=[descL '/mintemp.kdh']
    end
    tlapse=[descL timelapse];
    mkdir(tlapse);
    lDescriptors=rawLagDescriptorsLeft
    LagrangeDescriptors2movit
else     
    tags=['_' num2str(tinilag)]
    left=1
    descL=[desc 'C_' param '_w' num2str(0) tags  '-' num2str(real_pos) '-' num2str(reuse_links)];
    mkdir(descL);
    fileHeader=[descL '/header.kdh']
    tlapse=[descL timelapse];
    mkdir(tlapse);
    if vheader==2
        fileHmax=[descL '/maxtemp.kdh']
        fileHmin=[descL '/mintemp.kdh']
    end
    tlapse=[descL timelapse];
    mkdir(tlapse);
    lDescriptors=rawLagDescriptorsLeft
    LagrangeDescriptors2movit
    
end
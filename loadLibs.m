% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%Importing libraries
%Basic libraries with utils for image processing mainly
mpath='../matlab_code/';
addpath([mpath '/io/']);%my io functions
addpath([mpath '/NIFTI/']);%3rd party lib to handle nifti and analyze
addpath([mpath '/utils/']);%several utilities shared
addpath([mpath '/filter'])
addpath([mpath '/emgm'])

%Library to handle the EMB tracking format
addpath([ 'libEMB/utils'])
addpath(['libEMB'])
%Library to create Movit visualization files
addpath(['libMovit'])

%Library with mechanics related functions.
addpath(['libMechanics/'])

%Library for Stats
addpath(['libStats/'])

%Library for misc utils
addpath(['libUtils/'])
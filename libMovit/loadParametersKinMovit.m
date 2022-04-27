% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%We use the same parameters used in the calculation...
%loadParametersKinematics

%%%%%%% Parameters Normalization %%%%%%%%%%%%%%%%%%%
vheader=2;%version of the header 1-old (movit-817) 2-extraheaderslocal (movit-818) 3-new
mxper=95;
mnper=5;
%bins, actually it would be 32768 bins (positive part of int16, but i keep
%this for security -movit may not like all values-)
bins=32767;%0 -> 32767 positive...
indexNorm=0;%We normalized the values as indexes instead of physical magnitudes.
globalNorm=0;%if 1 the index of each step are normalized with the global max/min, but if 0 the indexes are normalized with the timestep max/min

temp_ave=Tave;
spat_ave=Xave;

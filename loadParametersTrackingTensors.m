% /* --------------------------------------------------------------------------------------
%  * File:    loadParametersTrackingTensors.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

%Low level configuration of the tracking to tensors processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PARAMETERS SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Structure taken from tracking to manage it as intermediate data.
redoDataStructure=0

%Averaging parameters
Tscale=20%2*sigma for temporal ave
Xscale=40%2*sigma for spatial ave
spatialAve=1;%flag to apply the spatial averaging (temporal averaging is always performed)

%Least Squares processing parameters
X=40%scale for samples for Least Square
maxNeigh=100%limits in the number of samples taken for Least Squares
minNeigh=5
onlySimilar=1;%Prunning outliers of samples for Least Squwares



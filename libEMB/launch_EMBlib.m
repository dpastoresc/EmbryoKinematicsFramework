% /* --------------------------------------------------------------------------------------
%  * File:    launch_EMBlib.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

%  Code seed imported from analysis of Benoit Lombardot (CNRS, ISC-PIF, France)

%Launching Eulerian calculations for EMB tracking
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if build_mov_F_E_R == 1
    
    'Generating tensors'
    generate_save_kinematics(datasetBioemergences, dataPath, Tscale,Xscale,Rscale,radialCorrection, spatialAve, X,maxNeigh,minNeigh, onlySimilar, tiniMat, tfinalMat)   
    disp(num2str(toc/60,'%.2f'));
end

disp([dataset ' : total time elapsed ' num2str(toc/60,'%.2f') ' min'])














% /* --------------------------------------------------------------------------------------
%  * File:    Tracking2vectorField_Tave.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

%  Code imported from analysis of Benoit Lombardot (CNRS, ISC-PIF, France)(CREA, Polytechnique) 

%Getting displacements field performing a gaussian-kernel averaging of the
%trajectories where tscale=2*sigma of the gaussian kernel

function [mov dmov]= Tracking2vectorField_Tave(cent,child,mother,t,Tave)
% cent: centers for each cell in t
% child: list of cell children identifiers for each cell in t
% mother: list of cell mother identifiers for each cell in t
% t:time step
% Tave: parameter os temporal averaging

%some dependent parameters that are fixed for any analysis.
%max move allowed
maxmoveallowed=9;%it was a magical number before
%min temporal integration
limitScale = Tave/3;%cells require a minimum of averaging to be considered 


ncell = size(mother{t},2);
cellList = 1:ncell;
tsigma = Tave/2; %Tave in the general documentation defiens the temporal averaging

%getting the chuncks of trajectories to use for filtering
%we go as far as [-Tave, Tave] and use Tave=2*sigma
tmin=max(1,t-Tave); tmax=min(length(mother),t+Tave+1);

%retrieve the trajectory chunck from the emb structures
[traj cellref] = get_trajs(cellList, mother, child, t, tmin, tmax);

%'checking'
%size(traj)
%storage of temporal averaged velocity field
mov = zeros(3,size(traj,2));
nmov = zeros(1,size(traj,2));

%storage of raw velocity field
rmov = zeros(3,size(traj,2));
nrmov = zeros(1,size(traj,2));

%perfoming filtering
for i=1:size(traj,1)-1
    %here we build the gaussian weights centered in on sample. This weights
    %will be used for each point of the trajectory (the interval around our
    %timepoint)
    tt = max(t-Tave,1)+i-1;     
    GWt = exp(-1*(t-tt).^2/(2*tsigma^2));%building the weigths
    
    %calculate the displacement
    %this corresponds to the online methods part of retrieving the vector
    %field
    dmov = cent{tt+1}(:,traj(i+1,:)) - cent{tt}(:,traj(i,:));
    dmovNorm = sum(dmov.^2,1).^0.5;%get the modulus

    sel = traj(i,:)>1 & traj(i+1,:)>1 & dmovNorm<maxmoveallowed;%here we constrain that we dont use samples so far away
    
    %calculating the average using the weights
    mov(:,sel) = mov(:,sel) + dmov(:,sel)*GWt;
    nmov(sel) = nmov(sel) + GWt;
end

%'checking'
%size(dmov)
%size(dmovNorm)
%size(mov)
%size(nmov)
%pause
%'ncell'
%ncell


%if there are new cells in this step (new dectections or mitosis)
if length(cellref)>ncell
    movaux = mov(:,ncell+1:end);
    nmovaux = nmov(:,ncell+1:end);
    cellrefaux = cellref(ncell+1:end);
    mov = mov(:,1:ncell);
    nmov = nmov(1:ncell);
    dmov = dmov(:,1:ncell);
    dmovNorm = dmovNorm(1:ncell);
    mov(:,cellrefaux) = mov(:,cellrefaux) + movaux;
    nmov(cellrefaux) = nmov(cellrefaux) + nmovaux;
end

%divigind by alpha=sum of weights.
mov = mov.*(ones(3,1)*(nmov>0))./max(ones(3,1)*nmov,1e-10);
mov = [mov;nmov]; 
mov(:,nmov<limitScale)= 0;%it alpha does not reach enough temporal steps
dmov = [dmov;dmovNorm];
%size(mov)
%size(nmov)
%size(dmov)
%size(dmovNorm)
%pause





% /* --------------------------------------------------------------------------------------
%  * File:    CreateVectorField_TensorField.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

%  Code imported from analysis of Benoit Lombardot (CNRS, ISC-PIF, France)(CREA, Polytechnique) 

% Filter trajectory field and get gradient of displacements tensor
% This is the core of the Eulerian decriptors


% NOTES:
% - speed is in microns/minute

function [Mov MeanMov TenseurGrad Neighbours Mov2order MeanMov2order RawMov] = Create_VectorField_TensorField(datas, t, params, my_flags)%radialCorrection, spatialAve)
%Mov vector field after getting vectors and temporal averaging
%MeanMov vector field after spatial averaging
%TenseurGrad F

t
'timestep within the array structure'
%S2(l) parameters... discretized function around the cell.
bins=10;
dl=4;
minSamples=3;

%Reading data.
cent = datas{1};%.cent;
mother = datas{2};%.child;
child = datas{3};%.mother;
spac = datas{4};%.spac;
sph = datas{5};%.sph;
tl = datas{6};%.tl;
tscale = params.tscale; % minutes = 2*Tsigma
Xscale = params.Xscale; % microns = 2*Xsigma
Rscale = params.Rscale; % microns
XmaxNeigh = params.XmaxNeigh;
XminNeigh = params.XminNeigh;
X=params.X;% en microns, this is the same of Xscale but for the derivatives calculation.
Xsigma = Xscale/2;
dt = tl(2)-tl(1); %time step in hours.
dtmin=dt*60; %in minutes
dtsec=dtmin*60; %in seconds

%Mov = zeros(size(cent{t}),'single');
%4values for meanMov
%23values for the tensor


size(cent)

%creating new data structures
MeanMov = zeros(4,size(cent{t},2),'single');
TenseurGrad = zeros(23,size(cent{t},2),'single');
%Neighbours = zeros(2,size(cent{t},2),'uint16');
%Mov2order = zeros(1,size(cent{t},2),'uint16');
%MeanMov2order = zeros(1,size(cent{t},2),'uint16');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% building displacement field + temp averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%David: Now we filter mitosis using a gaussian... change of cells and cell
%%mitosis are avoided!! (nice)
%%Mov is the mov of independent trajectories without considering neighbours
%%Averaging in time!!! we keep this averaging... 
Dt = ceil(tscale/(dtmin))
[Mov RawMov] = Tracking2vectorField_Tave(cent,child,mother,t,Dt); % moyenne des d�placement dans le temps
                                        % pond�ration gaussienne
                                        % si mitose moyenne des trajectoires ayant le m�me point de d�part
%'checking2'
%size(Mov)
%size(RMov)
%pause
                                       
%gaussian filtering weights (alpha in the doc). If 0, it has been discared...
Wt = Mov(4,:);
%displacements->microns/min!!! 
Mov = single(Mov(1:3,:)/dtmin);
RawMov = single(RawMov(1:3,:)/dtmin);

%if the entry was [0 0 0], the sample was discarded... we get the valid ones.
valid = sum(Mov.^2,1)>0;
svalid=size(valid);
sss=sum(valid(:));
discarded=svalid-sss
'valid eulerian samples, not integrated in the temporal averaging ones are discarded. im in movstraincalc2'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% S2(l) and continuum analysis for the raw velocity field%%%%%%%%%%%%%%%%%
%continuumDistribution_v2(Mov, cent{t}, t, 'mean', 0, 5, 10);
%pause
%NEW!!! We calculate the second order function                                        
[Mov2order, statsRaw]=secondOrderFunction_v3(Mov, cent{t}, dl, bins, minSamples, X);
%Mov2order=0
%'mov 2 order'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%spatial averaging%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% building neighbours
posList = cent{t};
neighList = [[0;0;0] , cent{t}(:,valid)];%some are discarded already
% size(posList)
% size(neighList)

%The radial correction uses sph to take less neighbours in the radian
%direction. Not enabled for now.
%if my_flags.radialCorrection==1
%    [neigh dist] = nn_getneigh3Sph(posList, neighList, Xscale, XmaxNeigh,spac,sph(1:3),Rscale);
%else
[neigh dist samplesAve] = buildNeighbours(posList, neighList, Xscale, XmaxNeigh,spac);
%end
neigh = neigh(2:end,:);%the first one is the number of neighs found

%Extra constrains for discarding divergent samples, pseudo-segmentation
dTheta_Max = pi/2;
minNormScale = 1/3;
minNorm=0.2;
% size(dist)
% max(dist(:))
% min(dist(:))
% pause

%creating weights
GWeight = exp(-1*dist.^2/(2*Xsigma^2));
GWeight(neigh==1)=0;
Movaux = [[0;0;0],Mov(:,valid)];
for i=2:size(Mov,2)
    mov0 = Mov(:,i);
    mov = Movaux(:,neigh(:,i));    
    %The slection of neighs with the extra contraints.
    sel = neigh_select(mov0,mov,dTheta_Max,minNormScale,minNorm);
    sel = sel' & neigh(:,i)>1;
    neigh_aux = neigh(sel,i);
    W = GWeight(sel,i);
    MeanMov(4,i) = sum(W.*Wt(neigh_aux)');
    
    %save neighbours for each cell. old
    %Neighbours(1,i)=size(find(sel>0),1);  
    if MeanMov(4,i)>0
        MeanMov(1:3,i) = sum(Movaux(:,neigh_aux).*(ones(3,1)*W'),2)./sum(W);
    end    
     
end

%%%%%%%%%% S2(l) and continuum analysis for the averaged velocity field%%%%%%%%%%%%%%%%%
[MeanMov2order, Neighbours]=secondOrderFunction_v3(MeanMov, cent{t},dl, bins, minSamples, X);
'mean mov 2 order'
%continuumDistribution_v1(MeanMov, cent{t}, t, 'mean', Xscale, dl, bins);
%pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%we use now X instead of Xscale obtained from the S2(l) function

%%David_ once we have averaged, we perform the calculation of derivatives.
%construction des d�formations moyennes aux �chelles Xscale, tscale
posList = cent{t};
neighList = cent{t};
ratio=Xscale/X%new ratio for max samples. review this
maxSamples=XmaxNeigh/ratio;
%Extra constrains for discarding neighbours. we have less strict conditions
%here as we want descriptors showing boundaries.
if my_flags.onlySimilar
  'same sample restrictions'
else
    'we take all samples for gradient'
    dTheta_Max = pi;
    minNormScale = 1/5;
    minNorm=0.1;
end
% dTheta_Max = pi/2; old ones.
% minNormScale = 1/3;
% minNorm=0.2;

%if my_flags.radialCorrection==1
%    [neigh dist] = nn_getneigh3Sph(posList, neighList, X, XmaxNeigh,spac,sph(1:3),Rscale);
%else
[neigh dist samplesGradient] = buildNeighbours(posList, neighList, X, maxSamples,spac);
%end
neigh = neigh(2:end,:);

%again weights, but they are not used.
GWeight = exp(-1*dist.^2/(2*X^2));
GWeight(neigh==1)=0;
%select if we want to have the averaged field or the raw field. we want of
%course the averaged field.
if my_flags.spatialAve==1
    FMovNorm = sum(MeanMov(1:3,:).^2,1).^0.5;
    FMov=MeanMov;
else
    FMovNorm = sum(Mov(1:3,:).^2,1).^0.5;  
    FMov=Mov;
end
for i=2:size(FMov,2)
    mov0 = FMov(1:3,i);
    mov = FMov(1:3,neigh(:,i));
    sel = neigh_select(mov0,mov,dTheta_Max,minNormScale,minNorm);
    sel = sel' & neigh(:,i)~=i & neigh(:,i)>1 & FMovNorm(neigh(:,i))'>0;
    %minimum number of samples
    if sum(sel)>XminNeigh
        neigh_aux = neigh(sel,i);
        size(neigh)
        size(sel)
        W = GWeight(sel,i);
        pos0 = cent{t}(:,i);
        pos  = cent{t}(:,neigh_aux);  
        
        dpos0 = FMov(1:3,i);
        dpos  = FMov(1:3,neigh_aux);       
       
        %Neighbours(2,i)=size(find(sel>0),1);
        %David: Here it does the operation with the derivatives!
        [F res resRel errCoeff]= get_gradientDeformationTensor(pos0, dpos0, pos, dpos, W, 'std');
        if size(F,2)==4
            shift = F(:,4);
        else
            shift = [0;0;0];
        end
        F = F(:,1:3);
        errCoeff = errCoeff(1:3,1:3);
        % F est le tenseur gradient: dx = F dX .
        % gradU = F-I est le gradient du champ de d�placement
        % E = 1/2 * ( ( F' * F ) - diag(ones(1,3))  );
        % E = 1/2 * ( gradU + gradU' + gradU'*grandU);
        % E est le tenseur des d�formations de Green-Lagrange
        TenseurGrad(:,i) = [ F(:) ; shift(:); res; resRel; errCoeff(:) ];
    end
end

%This function removes samples wich are not consistent with the sample
%reference.
%They are used to remove neighbours in the velocity field averaging and
%also in the gradient of displacements calculation
function sel = neigh_select(mov0,mov,dTheta_Max,NormScale,MinNorm)

n_neigh = size(mov,2);
norm0 = sum(mov0.^2).^0.5;
norm  = sum(mov.^2,1).^0.5;
if norm0>0
    selnorm = (norm0-norm)< max(MinNorm,(1-NormScale)*norm0) & (norm0-norm)> min(-MinNorm,(1-1/NormScale)*norm0);

    if norm0 > MinNorm
        cosTh = sum((mov0*ones(1,n_neigh)).*mov,1)./max(norm0.*norm,0.01);
        selTh = cosTh>cos(dTheta_Max);
    else
        selTh = ones(1,n_neigh)==1;
    end
    sel = selnorm & selTh;
else
    sel = ones(1,n_neigh)==1;
end


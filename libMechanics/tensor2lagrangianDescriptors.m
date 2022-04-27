function [Lag, LagLeft]=tensor2lagrangianDescriptors(globalTensor, T,countGoodSteps,countZeros, lagvelo)

% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%GETTING FINITE TIME DESCRIPTORS FROM TENSORS ALONG THE MATERIAL PARTICLES
%(LAGRANGIAN DESCRIPTORS)

%FTLE, FTLE isocoric, lagvelo, s, J, MIC1, MIC2, rot, rotVector, c1, ev1', c2, ev2', c3, ev3', shearAlphaAbs, shearAlphaIntermediate, countGoodSteps, countZeros
ndesc=27;
Lag=double(zeros(1,ndesc));
LagLeft=double(zeros(1,ndesc));
I=[1 0 0; 0 1 0; 0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating tensors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Forwards %%%%%%%%%%%%%%%%%%%
%Cauchy-green right tensor -> C
C=( globalTensor' * globalTensor );
%Cauchy-green left tensor -> B
B=( globalTensor * globalTensor' );
%Strain lagrange, strain accumulated seen from the original state
SL=0.5*(C-I);
%Strain euler, strain accumulated seen from the deformed state
SE=0.5*(I-inv(B));
%calculate U => F=RU
[R U]=F2RU(globalTensor);

%%%%%%%%%%%%%% Backwards %%%%%%%%%%%%%%%%%%
Fback=inv(globalTensor);
Cback=(Fback' * Fback);%???
Bback=(Fback * Fback');%???

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%J and Isocoric tensor. no volume gain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J=det(globalTensor);%=sqrt(det(C))

%%%% Fordwards
%Doing the isocoric from F can give us back complex eigen values
IF=(J^(-1/3))*globalTensor;
%Isochoric Cauchy-Green Tensors
IC=( IF' * IF );
IB=( IF * IF' );

%%%% Backwards
IFback=(det(Fback)^(-1/3))*Fback;
%D=(det(U)^(-1/3))*U;
ICback=( IFback' * IFback );
%Cauchy-green left tensor -> B
IBback=( IFback * IFback' );

%We dont use the SL and SE tensors, we just use the Cauchy-Green tensors
%Careful: no deformation -> Ci=1

%%%%%%%%%%%% EigenValues %%%%%%%%%%%%%%%%%%%%%
%%%% Strain in the origin
[CVectors, CValues]=eig(C);
CValues=[CValues(1,1) CValues(2,2) CValues(3,3)];
evaluesRank=orderEigValuesIndex_n1(CValues);
ce1=CValues(evaluesRank(1));
ce2=CValues(evaluesRank(2));
ce3=CValues(evaluesRank(3));
cv1=CVectors(:,evaluesRank(1));
cv2=CVectors(:,evaluesRank(2));
cv3=CVectors(:,evaluesRank(3));

%Isochoric version
[ICVectors, ICValues]=eig(IC);
ICValues=[ICValues(1,1) ICValues(2,2) ICValues(3,3)];
evaluesRank=orderEigValuesIndex_n1(ICValues);
ice1=ICValues(evaluesRank(1));
ice2=ICValues(evaluesRank(2));
ice3=ICValues(evaluesRank(3));

%Backwards version (normal and isochoric)
[CbVectors, CbValues]=eig(Cback);
CbValues=[CbValues(1,1) CbValues(2,2) CbValues(3,3)];
[ICbVectors, ICbValues]=eig(ICback);
ICbValues=[ICbValues(1,1) ICbValues(2,2) ICbValues(3,3)];
iceb1=ICbValues(evaluesRank(1));
iceb2=ICbValues(evaluesRank(2));
iceb3=ICbValues(evaluesRank(3));

%%%% Strain in the updated position
[BVectors, BValues]=eig(B);
BValues=[BValues(1,1) BValues(2,2) BValues(3,3)];
evaluesRank=orderEigValuesIndex_n1(BValues);
be1=BValues(evaluesRank(1));
be2=BValues(evaluesRank(2));
be3=BValues(evaluesRank(3));
bv1=BVectors(:,evaluesRank(1));
bv2=BVectors(:,evaluesRank(2));
bv3=BVectors(:,evaluesRank(3));

%Isochoric version
[IBVectors, IBValues]=eig(IB);
IBValues=[IBValues(1,1) IBValues(2,2) IBValues(3,3)];
evaluesRank=orderEigValuesIndex_n1(IBValues);
ibe1=IBValues(evaluesRank(1));
ibe2=IBValues(evaluesRank(2));
ibe3=IBValues(evaluesRank(3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vectors resolution to integer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
factor=30000;
bv1 = vector2longvector(bv1, factor);
bv2 = vector2longvector(bv2, factor);
bv3 = vector2longvector(bv3, factor);
cv1 = vector2longvector(cv1, factor);
cv2 = vector2longvector(cv2, factor);
cv3 = vector2longvector(cv3, factor);
%rotVector = vector2longvector(rotVector, factor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FTLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FTLE
leStretch=log(max(CValues(:)));
%Isocoric FTLE
leDis=log(max(ICValues(:)));
%FTLE backwards
leBStretch=log(max(CbValues(:)));
%Isocoric FTLE backwards
leBDis=log(max(ICbValues(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shear angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sIC=sort(ICValues(:));
shearAlphaAbs=(sqrt(sIC(3))-sqrt(sIC(1)))/2;
shearAlphaIntermediate=(sqrt(sIC(3))-sqrt(sIC(2)))/2;

% sICb=sort(ICbValues(:));
% shearAlphaAbsB=(sqrt(sICb(3))-sqrt(sICb(1)))/2
% shearAlphaIntermediateB=(sqrt(sICb(3))-sqrt(sICb(2)))/2
% cmm=max(CValues(:));
% cmi=min(CValues(:));
% shearA=(sqrt(cmm(1))-sqrt(cmi(1)))/2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shear angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MIC1=ice1+ice2+ice3;
MIC2=ice1*ice2+ice2*ice3+ice1*ice3;
% IC2=ce1*ce2+ce2*ce3+ce1*ce3
% MIC2B=iceb1*iceb2+iceb2*iceb3+iceb1*iceb3
% MIBC2=ibe1*ibe2+ibe2*ibe3+ibe1*ibe3
% T

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%theorem Euler. Euler-rodrigues
[rotVector rot]=getFiniteRotation(R);
rotVector = vector2longvector(rotVector, factor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=50;%Sign labelling start up
test=CValues<1;
s=s+sum(test(:))*50;

%%%% Intialize values %%%%%
Lag=[0 0 0 50 1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
LagLeft=[0 0 0 50 1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
FTLE=0;
FTLEiso=0;
veloAve=0;
FTLEB=0;
FTLEBiso=0;

if T>0%be carefull
    FTLE=leStretch/T;
    FTLEiso=leDis/T;
    FTLEB=leBStretch/T;
    FTLEBiso=leBDis/T;
    veloAve=lagvelo/T;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Saving descriptors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Careful with -inf values for  LCS (when value is zero)
    Lag(:,1)=FTLE;
    Lag(:,2)=FTLEiso;%sumai;%J-1;%leCompr;%/T;
    Lag(:,3)=veloAve;
    Lag(:,4)=s;
    Lag(:,5)=J;%or log(J)??
    Lag(:,6)=MIC1-3;
    Lag(:,7)=MIC2-3;
    Lag(:,8)=rot;
    Lag(:,9:11)=rotVector;
    Lag(:,12)=ce1;% or log... carefull this is not linear colormap
    Lag(:,13:15)=cv1';
    Lag(:,16)=ce2;
    Lag(:,17:19)=cv2';
    Lag(:,20)=ce3;
    Lag(:,21:23)=cv3';
    Lag(:,24)=shearAlphaAbs;
    Lag(:,25)=shearAlphaIntermediate;
    Lag(:,26)=countGoodSteps;
    Lag(:,27)=countZeros;
    
    %Left
    LagLeft(:,1)=FTLEB;
    LagLeft(:,2)=FTLEBiso;
    LagLeft(:,3)=veloAve;
    LagLeft(:,4)=s;
    LagLeft(:,5)=J;%or log(J)??
    LagLeft(:,6)=MIC1-3;
    LagLeft(:,7)=MIC2-3;
    LagLeft(:,8)=rot;
    LagLeft(:,9:11)=rotVector;
    LagLeft(:,12)=be1;
    LagLeft(:,13:15)=bv1';
    LagLeft(:,16)=be2;
    LagLeft(:,17:19)=bv2';
    LagLeft(:,20)=be3;
    LagLeft(:,21:23)=bv3';
    LagLeft(:,24)=shearAlphaAbs;
    LagLeft(:,25)=shearAlphaIntermediate;
    LagLeft(:,26)=countGoodSteps;
    LagLeft(:,27)=countZeros;
end
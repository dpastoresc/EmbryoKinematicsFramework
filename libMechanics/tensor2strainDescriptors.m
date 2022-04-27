function Descriptors=tensor2strainDescriptors(E)

% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%GETTING INCREMENTAL STRAIN DESCRIPTORS

%Parameters
%           - Tensor 3x3 E
%Out:
%           -  Descriptors
%               Qs
%               e1 / vn_e1: eigenvalues e1 strain
%               e2 / vn_e2: eigenvalues e2 strain
%               e3 / vn_e3: eigenvalues e3 strain
%               Qd
%               Qs-Qd
%               Max Shear Angle
%               de1 / vn_de1: eigenvalues e1 shear
%               de2 / vn_de2: eigenvalues e2 shear
%               de3 / vn_de3: eigenvalues e3 shear

%We dont want to store doubles now..
%We normalize the eigenvector and convert it to uint16.
vvaluenorm=30000;
n=28;
Descriptors=double(ones(1,n));
Descriptors=Descriptors.*(-1);

%Maybe the tensor is all zero. non valid
if ~any(E~=0)
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Eigenvalues and eigenvectors %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We have to look at this and see when the matrix is not symmetric... next improvement
[EigenVectors EigenValues]=eig(E);
%We order them by biggest to smalles (3,3 is the biggest)
ce=[EigenValues(1,1) EigenValues(2,2) EigenValues(3,3)];
%the biggest 1 -> to the smallest 3
evaluesRank=orderEigValuesIndex(ce);
se1=ce(evaluesRank(1));
se2=ce(evaluesRank(2));
se3=ce(evaluesRank(3));
v1=EigenVectors(:,evaluesRank(1));
v2=EigenVectors(:,evaluesRank(2));
v3=EigenVectors(:,evaluesRank(3));
e1=abs(se1);%we order it like this, so e1 is the biggest strain. We store the abs for colormap representation. We should keep sign too
e2=abs(se2);
e3=abs(se3);

%Qs
qs=e1*e2+e2*e3+e1*e3;
clear evaluesRank
clear EigenValues
clear EigenVectors
%Principal directions
nn=sqrt(v1(1)^2+v1(2)^2+v1(3)^2);
v1=v1./nn;
v1=int16(v1.*vvaluenorm);%This introduces some error in the directions...
nn=sqrt(v2(1)^2+v2(2)^2+v2(3)^2);
v2=v2./nn;
v2=int16(v2.*vvaluenorm);
nn=sqrt(v3(1)^2+v3(2)^2+v3(3)^2);
v3=v3./nn;
v3=int16(v3.*vvaluenorm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Deviatoric tensor to measure the rearrangement / shears
%(no expansion/compression)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equivalent to shear deformation
tracemean=(se1+se2+se3)/3;
T=[[1;0;0;] [0;1;0;] [0;0;1]].*tracemean;
%Deviatoric tensor: no expansion or compression
Dev=E-T;

[DEigenVectors DEigenValues]=eig(Dev);
cde=[DEigenValues(1,1) DEigenValues(2,2) DEigenValues(3,3)];
devaluesRank=orderEigValuesIndex(cde);
de1=cde(devaluesRank(1));
de2=cde(devaluesRank(2));
de3=cde(devaluesRank(3));
vd1=DEigenVectors(:,devaluesRank(1));
vd2=DEigenVectors(:,devaluesRank(2));
vd3=DEigenVectors(:,devaluesRank(3));

%Qd
qd=-(de1*de2+de2*de3+de1*de3);
clear devaluesRank
clear DEigenValues
clear DEigenVectors

%Qs-Qd
DistortionRatio=qs-qd;

%Maximum shear
MSA=(max(cde)-min(cde))/2;%(de1-de3)/2;

%Principal directions
nn=sqrt(vd1(1)^2+vd1(2)^2+vd1(3)^2);
vd1=vd1./nn;
vd1=int16(vd1.*vvaluenorm);
nn=sqrt(vd2(1)^2+vd2(2)^2+vd2(3)^2);
vd2=vd2./nn;
vd2=int16(vd2.*vvaluenorm);
nn=sqrt(vd3(1)^2+vd3(2)^2+vd3(3)^2);
vd3=vd3./nn;
vd3=int16(vd3.*vvaluenorm);

%to get it into the array
vd3=double(vd3');
vd2=double(vd2');
vd1=double(vd1');
v3=double(v3');
v2=double(v2');
v1=double(v1');

Descriptors=[qs se1 v1 se2 v2 se3 v3 qd DistortionRatio MSA de1 vd1 de2 vd2 de3 vd3];



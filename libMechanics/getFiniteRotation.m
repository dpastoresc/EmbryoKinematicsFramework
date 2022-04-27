% /* --------------------------------------------------------------------------------------
%  * File:    getFiniteRotation.m
%  * Date:    01/06/2016
%  * Author:  David Pastor Escuredo, david.pastor@upm.es
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015, David Pastor Escuredo - david.pastor@upm.es
%  All rights reserved.

function [axis angle]=getFiniteRotation(R)

%check R matrix is a valid rotation
%x=0.000001;
%R=[1 0 0; 0 cosd(x) -sind(x); 0 sind(x) cosd(x)];

det(R);
%should be zeros
%check=R'-inv(R);

%Euler-Rodrigues...
cosX=(trace(R)-1)/2;
angle=acosd(cosX);

%Eigenvalues and vectors to get the plane of rotation to get the axis.
[Rvectors Rvalues]=eig(R);

%check which eigenvector is the complex one
iv=abs(sum(imag(Rvalues)));
ivi=find(iv>0);
% iv
% pause

if numel(ivi)>0
    v1=Rvectors(:,ivi(1));
    %nn=sqrt(v1(1)^2+v1(2)^2+v1(3)^2);
    %v1=v1./nn; 
    v2=Rvectors(:,ivi(2));
    %nn2=sqrt(v2(1)^2+v2(2)^2+v2(3)^2);
    %v2=v2./nn2; 
    axis=cross(v1, v2);
    nn=sqrt(axis(1)^2+axis(2)^2+axis(3)^2);
    axis=axis/nn;
    axis=axis';

else    
   axis=[0 0 0];
   angle=0;
end

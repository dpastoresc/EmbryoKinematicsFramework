function Descriptors=tensor2topologyDescriptor(F, velo, mvelo, neighs)

% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%GETTING FLOW TOPOLOGY DESCRIPTORS
%           - Tensor H: 3x3
%Out:
%           - Descriptors=[speed speedAve mvelo q top p D e exp_nv c comp_nv angle rota_nv]

%%%% Init to -1 (control error value)
vvaluenorm=30000;
Descriptors=double(ones(1,21));
Descriptors=Descriptors.*(-1);

%%%% Velocity
speed=double(sqrt(velo(1)^2+velo(2)^2+velo(3)^2));
speedAve=double(sqrt(mvelo(1)^2+mvelo(2)^2+mvelo(3)^2));
mvelo=mvelo(1:3)./speedAve;
mvelo=int16(mvelo.*vvaluenorm);
mvelo=double(mvelo');
neighs=double(neighs);

%Maybe the tensor is all zero->not enough samples, error in the Least
%Squares
%Non valid->then all is -1
if ~any(F~=0)
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRADIENT OF DISPLACEMENTS %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We operate with H=F-I
I=[1 0 0; 0 1 0; 0 0 1];
H=F-I;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Eigenvalues and eigenvectors %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We have to look at this and see when the matrix is not symmetric... next improvement
[EigenVectors EigenValues]=eig(H);
e1=EigenValues(1,1);
e2=EigenValues(2,2);
e3=EigenValues(3,3);
v1=EigenVectors(:,1);
v2=EigenVectors(:,2);
v3=EigenVectors(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Analysis of linear components of topology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% P
%P>0 expands                %P<0 compresses
%p=real(e1)+real(e2)+real(e3);%P invariant
p=trace(H);

%%%% Topology
s=50;%Sign labelling start up

%%%% Planor compression or expansion P+Topology
e=0;
c=0;
%Number of eigenvectors stored
vcountn=1;%Vector for plane of compresion
vcount=1;%Vector of plane of expansion
%check eigenvalues
if real(e1)<0
    s=s+50;%anytime we find a compressive eigenvalue we change the label
    c=c+real(e1);
    planen(:,vcountn)=v1;
    vcountn=vcountn+1;
else
    e=e+real(e1); %if expansion
    plane(:,vcount)=v1;
    vcount=vcount+1;
end
if real(e2)<0
    c=c+real(e2);
    s=s+50;
    planen(:,vcountn)=v2;
    vcountn=vcountn+1;
else
    e=e+real(e2);
    plane(:,vcount)=v2;
    vcount=vcount+1;
end
if real(e3)<0
    c=c+real(e3);
    s=s+50;
    planen(:,vcountn)=v3;
    vcountn=vcountn+1;
else
    e=e+real(e3);
    plane(:,vcount)=v3;
    vcount=vcount+1;
end

%Compressing plane angle
exp_nv=[0;0;0];
if(s==100)
    %Calculate the vector of the surface of expansion
    exp_nv=cross(plane(:,1),plane(:,2));
    nn=sqrt(exp_nv(1)^2+exp_nv(2)^2+exp_nv(3)^2);
    exp_nv=exp_nv./nn;
    exp_nv=int16(exp_nv.*vvaluenorm);
else
    e=-1;
end
exp_nv=double(exp_nv');

%Compressing plane angle
comp_nv=[0;0;0];
if(s==150)
    c=c*-1;%magnitud is positive!
    %Calculate the vector of the surface of compression
    comp_nv=cross(planen(:,1),planen(:,2));
    nn=sqrt(comp_nv(1)^2+comp_nv(2)^2+comp_nv(3)^2);
    comp_nv=comp_nv./nn;
    comp_nv=int16(comp_nv.*vvaluenorm);
else
    c=-1;
end
comp_nv=double(comp_nv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 i1=imag(e1);
%                 i2=imag(e2);
%                 i3=imag(e3);
%                 icount=1;
%                 if(abs(i1)>0)
%                   planeRotation(:,icount)=v1;  icount=icount+1; end
%                 if(abs(i2)>0)
%                   planeRotation(:,icount)=v2;  icount=icount+1;   end
%                 if(abs(i3)>0)
%                   planeRotation(:,icount)=v3;  icount=icount+1; end
%
%                 %planeRotation
%                 rota_nv=[0;0;0];
%                 if(icount>1)
%                     %Calculate the vector of the surface of compression
%                     rota_nv=cross(planeRotation(:,1),planeRotation(:,2));
%
%                     nn=sqrt(rota_nv(1)^2+rota_nv(2)^2+rota_nv(3)^2);
%                     rota_nv=rota_nv/nn;
%                     rota_nv=int16(rota_nv*vvaluenorm);
%                     %nn=sqrt(rota_nv(1)^2+rota_nv(2)^2+rota_nv(3)^2)
%                     rota_nv=double(rota_nv');
%                     %Rotation ratio
%                     rat1=abs(imag(e1)/(real(e1)+0.000000000001));
%                     rat2=abs(imag(e2)/(real(e2)+0.000000000001));
%                     rat3=abs(imag(e3)/(real(e3)+0.000000000001));
%                     rat=max([rat1 rat2 rat3]);
%                     if rat<1
%                         rat=0;%If there is no more rotation than linear movement
%                     end
%                     %Rotation component
%                     rota=max([abs(imag(e1)) abs(imag(e2)) abs(imag(e3))]);
%                 else
%                     rat=-1;
%                     rota=-1;
%                 end
%               es=[e1 e2 e3]
%               me=max(abs(es));
%               ime=find(abs(es)==me);
%               em=es(ime);
%Alternative way to calculate the rotation. Should be
%compared...
%W = (H-H')/2;           %Spin tensor
%w1=W(3,2);
%w2=W(1,3);
%w3=W(2,1);
%wv=[w1 w2 w3];
%wspeed=sqrt(w1^2+w2^2+w3^2);
%we apply the general decomposition of R to get the
%rotation... we could use the linear one for the eulerian
%but it does not matter...
[R U] = F2RU(F);
%                 R
%                 EigenVectors
%                 EigenValues
[axis angle]=getFiniteRotation(R);

rota_nv = double(vector2longvector(axis,vvaluenorm));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=1/2*((trace(H)^2)-trace(H*H));
q=e1*e2+e1*e3+e2*e3;
% W=1/2*(H-H')
% trace(W)
% [EigenVectors EigenValues]=eig(W);
% x1=EigenValues(1,1)
% x2=EigenValues(2,2)
% x3=EigenValues(3,3)
% x1+x2+x3
% e1
% e2
% e3
% qprima=x1*x2+x1*x3+x2*x3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation Discriminant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=-det(H);
p1=-trace(H);
D=27*r^2+(4*p1^3 - 18*p1*q)*r+(4*q^3-p1^2*q^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% topology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
top=0;
%50 expands and rotates
%100 expands and no rotates
%150 compress and rotates
%200 compress and no rotates
if p>0 
    if D>0
        top=50;%values for a colormap 255.
    else
        top=100;
    end
else
    if D>0
        top=150;
    else
        top=200;
    end
end

s=double(s);
Descriptors=double([speed speedAve mvelo top p q D e exp_nv c comp_nv angle rota_nv]);% rat]);



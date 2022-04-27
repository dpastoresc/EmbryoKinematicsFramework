% /* --------------------------------------------------------------------------------------
%  * File:    vf2migrationDescriptors.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

function Descriptors=vf2migrationDescriptors(velo, mvelo, rvelo, velosl, mvelosl, neighs, oldneigh, cellDiam, convergence)

% Descriptors field
% Descriptors=[speed velo(3) speedAve veloAve(3) speedR veloR(3) 2orderMovMod 2orderMovMeanComplete(3) neighXave neighGainXave cellDiam Convergence]

            %We dont want to store doubles now..
            %We normalize the eigenvector and convert it to uint16.
            vvaluenorm=30000;
            n=7;
            %Descriptors=double(zeros(1,n));
            %escriptors=Descriptors.*(-1);
         
            speed=double(sqrt(velo(1)^2+velo(2)^2+velo(3)^2));
            speedAve=double(sqrt(mvelo(1)^2+mvelo(2)^2+mvelo(3)^2));   
            speedR=double(sqrt(rvelo(1)^2+rvelo(2)^2+rvelo(3)^2));   
            
            %velosl
            %mvelosl
            %neighs
            neighGain=neighs-oldneigh;
            %pause
            
            velo=velo(1:3)./speed;
            velo=int16(velo.*vvaluenorm);  
            velo=double(velo');
          
            mvelo=mvelo(1:3)./speedAve;
            mvelo=int16(mvelo.*vvaluenorm);  
            mvelo=double(mvelo');
            
            rvelo=rvelo(1:3)./speedR;
            rvelo=int16(rvelo.*vvaluenorm);  
            rvelo=double(rvelo');
           
            neighGain=double(neighs-oldneigh);          
            neighs=double(neighs);
            Descriptors=[speed velo speedAve mvelo speedR rvelo velosl(4) double(mvelosl') neighs(1) neighGain(1) cellDiam convergence];
            Descriptors(Descriptors==-1)=-100;%we have this descarding value here...
    

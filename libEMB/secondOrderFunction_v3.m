% /* --------------------------------------------------------------------------------------
%  * File:    secondOrderFunction_v3.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

%v2 fixes a lot of stuff and exports the results for .mat

%STATS in big table / ring
%1-meanVel 2-meanVel_dif 3- stdVel 4-meanAngleRef 5-stdAngleRef_dif 6-stdAngleRef 
%7-disMean 8-disStd 9-neigh 10-density 11-meanAngle 12- stdAngle 
%13- pow (all rings the same) %14-C2 (all rings the same)

function [SecOrderFunc, SpatialStats]=secondOrderFunction_v3(Mov, posList, dl, nbins, minSamples, X)

     
 %Mov={u,v,w} for all cells (postList and Mov have coherent indexes)
 %postList={x,y,z} for all cells
 %dl is the length difference to measure the convergence
 %nbins the number of intervals we use to approximate the decay
 %minsamples is the min number of cells we need to use one of the bins
  if nargin < 3
    dl = 6
 end  
 if nargin < 3
     dl = 6
 end
 if nargin < 4
    nbins = 10
 end
 if nargin < 5
    minSamples=3
 end
  if nargin<6
     X=20;
 end
 
 l=[dl*nbins:-dl:dl]; %upper limits of the rings
 ncells=size(posList,2);%number of cells in this timestep
 SecOrderFunc=zeros(4,ncells);%the descriptor is the powerlaw (m,b) y=b*x^m
 SpatialStats=zeros(5,ncells); 
 
 l
 %nbins
 Xring=floor(X/dl);
 Xring=nbins-Xring+1;%inverted
 
for i=1:ncells   
    
    %cell ref values
    x=posList(1,i);
    y=posList(2,i);
    z=posList(3,i); 
    u=Mov(1,i);
    v=Mov(2,i);
    w=Mov(3,i);
    %cell ref speed
    m=sqrt(u^2+v^2+w^2);
    nv=[u v w]./m;
 
    %distances between ref and other cells
    dx=(posList(1,:)-x).^2;
    dy=(posList(2,:)-y).^2;
    dz=(posList(3,:)-z).^2;
    dAll=sqrt(dx+dy+dz);
    %mean(dAll)%     min(dAll)%     max(dAll)%     std(dAll)
    
    pdecayall=zeros(4,length(l));
    samples=zeros(1,length(l));%number of cells found within the n*dl
    anglesrelall=-1;
    daveall=-1;
    numOfSamples=0;
    
    c=1;%counter    
    for nl=1:length(l)
        
        %stats to save
        anglesrel=zeros(0);
        acount=1;
        dave=zeros(0);
        dcount=1;
        %%%%ring spatial radius (microns)
        l_up=l(nl);
        l_low=l_up-dl;
        %find the points at (n-1)*dl<d<n*dl in our centers list
        %this can be nice to measure tissue density...
        
        %To Do: improve window!!!
        isel=find(dAll>l_low & dAll<=l_up);
        %%%we keep the original selection=number of neighbours. We may
        %%%prune this selection statistically
        iselkeep=isel;
        
        %go through the array of cells in each cell
        if length(isel>0)
            mAll=sqrt(Mov(1,isel).^2+Mov(2,isel).^2+Mov(3,isel).^2);            
            NM=Mov(1:3, isel);
            selPos=posList(:,isel);
          
            for ir=1:length(isel)                
               
                %each sample inside the ring
                nvv=NM(:,ir)./norm(NM(:,ir));
                
                %checking the closes neighbour
                s1=selPos(1,:)-selPos(1,ir);
                s2=selPos(2,:)-selPos(2,ir);
                s3=selPos(3,:)-selPos(3,ir);
                dd=sqrt(s1.^2+s2.^2+s3.^2);
                idd=find(dd>0);
                dd=dd(idd);
                if length(dd)>0
                    [mdd imdd]=min(dd);
                    dave(dcount)=mdd;
                    dcount=dcount+1;
                    
                    %angle to the closes neigh
                    nvc=NM(:,imdd)./norm(NM(:,imdd));
                    ac=acosd(dot(nvv', nvc));
                    %nan values are omitted
                    if ~isnan(ac)
                       anglesrel(acount)=ac;
                       acount=acount+1;
                    end
                end                        
                
            end
            
        if nl>=Xring
            anglesrelall=horzcat(anglesrelall, anglesrel);
            daveall=horzcat(daveall, dave);
        end
        else          
            NM=Mov(1:3,:);
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Relative descriptors to the cell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %S2(l): second order structure function
        du=(NM(1,:)-u).^2;
        dv=(NM(2,:)-v).^2;
        dw=(NM(3,:)-w).^2;
        dvel=sqrt(du+dv+dw);
        dvel=dvel.^2;%square of the modulus of the vectorial difference.
        dall=horzcat(du', dv', dw', dvel');
        
        %pair distribution function        
        %dvel=abs(mAll-m);
        n= length(iselkeep);%we may have pruned the samples speed, but the number would  hold    
        numOfSamples=numOfSamples+n;
            
        if n>0
            ave=sum(dall,1)/n;%sum all the rows by column
            pdecayall(:,c)=ave';
            samples(c)=n;
        else
            samples(c)=0;
            pdecayall(:,c)=[0 0 0 0]';
            %we remove it later.
            %"we exclude d=0, so we can found no cells between in [0 dl]
        end
        c=c+1;       
    end       
    %pdecayall
    %l
    %samples
    %if there were no detections within n*dl we have to remove the bin to
    %avoid singularity
    %exceptValue. Propagate....
    powlaw=[-100 -100];
   % size(samples)
    for ic=1:4
       
        pdecay=pdecayall(ic,:);
        lx=l(Xring:nbins);
        samplesx=samples(Xring:nbins);
        pdecay=pdecay(Xring:nbins);
        
        %no cells found to make the analysis
        while pdecay(end)==0 & length(pdecay)>1                          
           pdecay=pdecay(1:end-1);
           lx=lx(1:end-1);
           samplesx=samplesx(1:end-1);           
        end
        %at least minSamples to make it more robust (one or 2 samples could
        %diverge for a false detection)
        while samplesx(end)<minSamples & length(pdecay)>1  
           pdecay=pdecay(1:end-1);
           lx=lx(1:end-1);
           samplesx=samplesx(1:end-1);
        end
        if length(pdecay)<3
            %'short time'
            
            SecOrderFunc(ic,i)=powlaw(1); 
            %is a singularity of this... shouldnt exist a powerlaw with 0 
            %if we dont have enough values to calculate we use 0        
        else
            %'getting the powlaw fit'
            %powlaw(1) should be 2 if the decay is governed l^2 
            %e^powlaw(2) is the multiplicative constant
            lx=lx-(dl/2);%we use the half of the bins division instead of the extrems
            powlaw = polyfit(log(lx), log(pdecay), 1);
            SecOrderFunc(ic,i)=powlaw(1);%we keep only the power exponent fore each component
            %powlaw
        end
    end
    %we save some stats
    %mean, std distance
    %mean, std of angles inside X
    %number of cells inside X
    
    %angles to the closest neighbour
    %anglesrelall
    %daveall
    [anglesrelall iar]=m_removeOutliers(anglesrelall, 90, 0, -1);
    [daveall idave]=m_removeOutliers(daveall, 98, 2, -1);    
    %anglesrelall
    %daveall
    md=mean(daveall);
    sd=std(daveall);
    ma=mean(anglesrelall(:));
    sa=std(anglesrelall(:));
    if md==-1
        sd=-1;
    end
    if ma==-1
        sa=-1;
    end
    desc=[md sd ma sa numOfSamples];
    SpatialStats(:,i)=desc;
    %pause
        
    
end

%['BigTable' tag  '_t' num2str(t) '_dl' num2str(dl) '_bins' num2str(nbins)] 
%pathh='../MechanicsData/stats/';
%save([pathh 'BigTable' tag '_t' num2str(t) '_dl' num2str(dl) '_bins' num2str(nbins) '_Xave' num2str(Xave)], 'BigTable');
%save([pathh 'S2l' tag '_t' num2str(t) '_dl' num2str(dl) '_bins' num2str(nbins) '_Xave' num2str(Xave)], 'accpdecay');
%save([pathh 'cells_t' num2str(t) '_dl' num2str(dl) '_bins' num2str(nbins) '_Xave' num2str(Xave)], 'cellstats');
%VisualizeStatsContinuum_v2(BigTable)
        
    
    
    
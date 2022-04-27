%David Pastor Escuredo (BIT-UPM / MAE-UCSD)
%2014
%Mechanics Framework
%(C) All rights reserved
%Having the velocity field (averaged or not) and the centers for each
%timestep, we calculate the second order function

%v2 outputs the powlaw for each of the components of the velocity field
function [SecOrderFunc samples l] = secondOrderFunction_v2(Mov, posList, dl, nbins, minSamples)

     
 %Mov={u,v,w} for all cells (postList and Mov have coherent indexes)
 %postList={x,y,z} for all cells
 %dl is the length difference to measure the convergence
 %nbins the number of intervals we use to approximate the decay
 %minsamples is the min number of cells we need to use one of the bins
 if nargin < 5
    minSamples=4
 end
 if nargin < 4
    nbins = 10
 end
 if nargin < 3
    dl = 5
 end  
 
 ncells=size(posList,2);%number of cells in this timestep
 SecOrderFunc=zeros(4,ncells);%the descriptor is the powerlaw (m,b) y=b*x^m
 %we keep m for x, y, z components and the module.
 
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
    
    %distances between ref and other cells
    dx=(posList(1,:)-x).^2;
    dy=(posList(2,:)-y).^2;
    dz=(posList(3,:)-z).^2;
    dAll=sqrt(dx+dy+dz);
    %size(dAll)
    %min(dAll)
    %max(dAll)

    %bins for the power law.
    %dl=5microns, should adjust it somehow with the min dist maybe
    %l->0 
    %dl=5
    %nbins=10
    l=[dl*nbins:-dl:dl];
    pdecay=zeros(4,length(l));%metric that should fit the powerlaw, 3 comp and mod
    samples=zeros(size(l));%number of cells found within the n*dl

    c=1;%counter    
    for nl=1:length(l)
        l_up=l(nl);
        l_low=l_up-dl;
        %find the points at (n-1)*dl<d<n*dl in our centers list
        %this can be nice to measure tissue density...
        isel=find(dAll>l_low & dAll<=l_up);
        %isel2=find(dAll<l_up);
    
        %speed of all cells... maybe we should use components better 
        %or even use an angle criterium?
        mAll=sqrt(Mov(1,isel).^2+Mov(2,isel).^2+Mov(3,isel).^2);
        %meanVel=mean(mAll(:))
        %minVel=min(mAll(:))
        %maxVel=max(mAll(:))
        
        %Mov is the speed vector for all the centers
        %size(mAll)
        %this if we want to make it by components.
        du=(Mov(1,isel)-u).^2;
        dv=(Mov(2,isel)-v).^2;
        dw=(Mov(3,isel)-w).^2;
        dvel=sqrt(du+dv+dw);
        dvel=dvel.^2;%square differences of velocity modulus??
            
        dall=horzcat(du', dv', dw', dvel');
        %dvel=abs(mAll-m);
        n=size(dvel,2);   
        if n>0
            ave=sum(dall,1)/n;%sum all the rows by column
            %size(dall)
            %size(ave)
            %dall
            %pause
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
    for ic=1:4
        pdecay=pdecayall(ic,:);
        lx=l;
        samplesx=samples;
        upd=unique(pdecay);
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
            SecOrderFunc(ic,i)=0; 
            %is a singularity of this... shouldnt exist a powerlaw with 0 
            %if we dont have enough values to calculate we use 0
        
        else
            %powlaw(1) should be 2 if the decay is governed l^2 
            %e^powlaw(2) is the multiplicative constant
            lx=lx-(dl/2);%we use the half of the bins division instead of the extrems
            powlaw = polyfit(log(lx), log(pdecay), 1);
            SecOrderFunc(ic,i)=powlaw(1);%we keep only the power exponent fore each component
        end
    end
end

    
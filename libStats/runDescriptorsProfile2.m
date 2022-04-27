% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%%%%%%%%%%%%%%%%% STATS ALONG TIME INTERVAL %%%%%%%%%%%%%%%%%%%%%%%%


%init data
%metric=zeros(t_limit-t_start+1, length(selections),13);
metric=[]

for tcurrent=t_start:t_limit
    
    'Processing step'
    tcurrent  
   
    %Obtains a profile for each selection (see loadStatsConfig.m)    
    for si=1:length(selections)
        %%%% Select time frame data
        sDesc=Desc(Desc(:,3)==tcurrent,:);
        snDesc=Desc(Desc(:,3)==(tcurrent+step),:);%if we want to get derivatives of descriptors in space
        sizeD=size(sDesc);
        
        %%%% Remove bad stuff
        sDesc=sDesc(sDesc(:,7)>invalidValue,:);
        snDesc=snDesc(snDesc(:,7)>invalidValue,:);
        sizeD_removed=size(sDesc);

        
        %%%% Process each selection
        ss=selections(si)
        %if purely material you would not need this, just the ids
        if spatial_mask
            select_t=selec(selec(:,3)==tcurrent,:);
        else
            %MAY NOT WORK
            %make material propagating from a tref
            if tref_propagate>-1
                select_t=selec(selec(:,3)==tref,:);
            else
                select_t=selec;
            end
        end
        all_selections_step=size(select_t); 
        ids=unique(select_t(select_t(:,2)==ss,1));
        
        %%%% Take the samples inside the domain selected
        ssDesc=sDesc;
        selec2= arrayfun(@(x) find(sDesc(:,1)==x), ids, 'UniformOutput', false);
        sids=size(selec2);
        selec2= cell2mat(selec2);
        sDesc=sDesc(selec2,:);
        entries=size(sDesc(:,1),1);
        entries_dif=size(unique(sDesc(:,1)),1);
        
        if debug_this
            ix=find(sDesc(:,descriptor)==-1);
            if size(ix,1)>0
                'ERROR: found -1'
                ix
                selec2(ix)
                ssDesc(selec2(ix),6:9)
                sDesc(ix,6:9)
                pause
            end
        end
        clear ssDesc
                   
        %%%%% Take Descriptor Selection
        Im=sDesc(:,descriptor);     
       
        %%%%% If gradient of spatial info->do it as Image
        if derive
            %%%% Descriptors to image
            Im=zeros(sizePixel);
            mask1=zeros(sizePixel);

            for i=1:size(sDesc,1)%first is control      
                xyz=uint16(sDesc(i,4:6)./spacingPixel(1:3))+1;
                xyz=min(uint16(xyz), uint16(sizePixel(1:3)));
                Im(xyz(1),xyz(2),xyz(3))=sDesc(i,descriptor);
                mask1(xyz(1),xyz(2),xyz(3))=1;
            end
            pix1=sum(mask1(:));%samples in the image

            %The derivatime in time of spatial information is done through
            %images

            Imn=zeros(sizePixel);
            mask2=zeros(sizePixel);
            for i=1:size(snDesc,1)
                xyz=uint16(snDesc(i,4:6)./spacingPixel(1:3))+1;
                %spatial tolerance in sparse data
                %ToDo: Should interpolate
                xyzm=xyz-spatialRad;
                xyzM=xyz+spatialRad;
                xyzM=min(xyzM, sizePixel(1:3));
                xyzm=max(xyzm, sizePixel(1:3));
                Imn(xyzm(1):xyzM(1),xyzm(2):xyzM(2),xyzm(3):xyzM(3))=snDesc(i,descriptor);
                mask2(xyzm(1):xyzM(1),xyzm(2):xyzM(2),xyzm(3):xyzM(3))=1;
            end
            mask1=mask1.*mask2;
            clear mask2;
            pix2=sum(mask1(:))
            Im=Im.*mask1;%make coincidence (compensate propagation)
            Imn=Imn.*mask1;
            Im=(Im-Imn).^2;
            Im=Im(mask1>0);        
            if debug_thus
                pix1
                pix2               
            end
            clear mask1;
            clear Imn;
        end
        
        %Data Before statistical cropping
        samplesTaken=length(Im);  
        if debug_this
            ixxx=find(Im==-1);
            if size(ixxx,1)>0
                'ERROR value -1'
                pause
            end            
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Remove outliars   
        
        %Check hist
        if make_hist
            hist(Im)
            title(['Histogram selection ' num2str(si) ' step ' num2str(tcurrent)])
            %save
            pause
        end
        
        mx=prctile(Im, pmax(DescriptorIndex));
        mn=prctile(Im, pmin(DescriptorIndex));
        if removeOutliars
            Im=Im(Im>=mn);
            Im=Im(Im<=mx);
        else
            Im(Im>mx)=mx;
            Im(Im<mn)=mn;
        end
        samplesInside=length(Im);
        
        %Check hist after
        if make_hist
            figure
            title(['Histogram (after outliers) selection ' num2str(si) ' step ' num2str(tcurrent)])
            %save
            pause
        end     
        
        if debug_this
            sids
            entries
            entries_dif
            sizeImP=size(Im)
            samplesTaken
            samplesInside
            pause
        end        
        metric(tcurrent+1,1,si)=mean(Im(:));
        metric(tcurrent+1,2,si)=median(Im(:));
        metric(tcurrent+1,3,si)=std(Im(:));
        metric(tcurrent+1,4,si)=mx;%max and min
        metric(tcurrent+1,5,si)=mn;
        metric(tcurrent+1,6,si)=prctile(Im, 99);
        metric(tcurrent+1,7,si)=prctile(Im, 1);
        metric(tcurrent+1,8,si)=prctile(Im, 75);
        metric(tcurrent+1,9,si)=prctile(Im, 25);
        metric(tcurrent+1,10,si)=entries;%number of entries
        metric(tcurrent+1,11,si)=entries_dif;%number of different entries check
        metric(tcurrent+1,12,si)=samplesTaken;%data taken
        metric(tcurrent+1,13,si)=samplesInside;%data after outliars
        metric(tcurrent+1,14,si)=ss;%the actual selection number
        metric(tcurrent+1,15,si)=tcurrent;%the actual selection number
     
    end
end
save(stats_file, 'metric');




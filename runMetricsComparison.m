% /* --------------------------------------------------------------------------------------
%  * File:    runMetricsComparison.m
%  * Date:    01/06/2015
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2013-2017, David Pastor Escuredo
%  with Biomedical Image Technology, UPM (BIT-UPM)
%  with BioEmergences, CNRS
%  with LifeD lab
%  All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% SPATIO-TEMPORAL CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%
addpath('colormaps')
cs=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Calculate distance to LBPs of the reference %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%D=zeros(size(clus,2), size(clus,2));
%Df=zeros(size(clus,2), size(clus,2));
D=zeros(size(clus,2)-5, 5);
Df=zeros(size(clus,2)-5, 5);

%WE COMPARE AT LEAST 3 TIME STEPS
%We skip t=1 as the lagrangian is all 1

for t=4:size(clus,2)-2
    % for t2=1:size(lbpRef,2)
    %min of distance to ref lbps
    lbpt=lbpRef(:,2:t);
    clust=clus(:,2:t);
    
    lbpf=lbpRef(:,t:end);
    clusf=clus(:,t:end);
    
    
    for j=1:size(clus,1)
        clustj=clust(j,:);
        clusfj=clusf(j,:);
        
        for b=1:size(lbpt,1)
            dclustj(b)=pdist2(clustj,lbpt(b,:),distanceT);
            dclusfj(b)=pdist2(clusfj,lbpf(b,:),distanceT);
        end
        dclust(j)=min(dclustj);
        dclusf(j)=min(dclusfj);
    end
    dref=prctile(dclust,25);
    dreff=prctile(dclusf,25);
    dref2=prctile(dclust,75);
    dreff2=prctile(dclusf,75);
    drefmin=min(dclust);
    dreffmin=min(dclusf);
    drefmax=max(dclust);
    dreffmax=max(dclusf);
    drefm=mean(dclust);
    dreffm=mean(dclusf);
    
    D(t-3,:)=[drefm dref dref2 drefmin drefmax];
    Df(t-3,:)=[dreffm dreff dreff2 dreffmin dreffmax];
    % end
end

dlmwrite([desc 'past_score' '-' tagDesc '.csv'],D)
dlmwrite([desc 'future_score' '-' tagDesc '.csv'],Df)

% maxsD=[0 0 0 0 0.0002 0.006 0.006 0.007 0]
% max(D(:))
% max(Df(:))
% 
% h=figure
% colormap('bone')
% %Df=Df./maxsD(DescriptorIndex);
% 
% if normType==1
%     
%     D=D./maxsD(DescriptorIndex);
%     D=1-D;
%     tagn='1'
% else
%     'hols'
% end
% 
% imagesc(D, [0 1])
% axis off
% screen_size = get(0, 'ScreenSize');
% set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
% if savePlots
%     saveas(h,[StatsPath dataset tagfolder filesep 'past_score' tagn '-' datref '-' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
% end
% 
% h2=figure
% colormap('pink')
% %Df=max(Df(:))-Df;
% 
% if normType==1
%     Df=Df./maxsD(DescriptorIndex);
%     Df=1-Df;
%     tagn=''
% else
%     'hola'
% end
% imagesc(Df, [0 1])
% axis off
% screen_size = get(0, 'ScreenSize');
% set(h2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
% if savePlots
%     saveas(h2,[StatsPath dataset tagfolder filesep 'future_score' tagn '-' datref '-' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
% end
% 
% max(D(:))
% max(Df(:))
% 
% pause
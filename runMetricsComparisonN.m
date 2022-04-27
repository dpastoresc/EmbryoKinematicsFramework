% /* --------------------------------------------------------------------------------------
%  * File:    runMetricsComparisonN.m
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
cs=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Calculate distance to LBPs of the reference %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%D=zeros(size(clus,2), size(clus,2));
%Df=zeros(size(clus,2), size(clus,2));
D=[]
Df=[]
%for j=1:size(clus,1)
for t=3:size(clus,2)-1
    for t2=3:size(lbpRef,2)-1
        %min of distance to ref lbps
        lbpt=lbpRef(:,2:t2);
        clust=clus(:,2:t);
        
        lbpf=lbpRef(:,t2:end);
        clusf=clus(:,t:end);
        
        
        for j=1:size(clus,1)
            clustj=clust(j,:);
            clusfj=clusf(j,:);
            
            lbpt_2=resizem(lbpt,[size(lbpt,1) size(clustj,2)],'bilinear');
            lbpf_2=resizem(lbpt,[size(lbpf,1) size(clusfj,2)],'bilinear');
            for b=1:size(lbpt,1)
                xt=squeeze(lbpt_2(b,:));
                dclustj(b)=pdist2(clustj,xt,distanceT);
            end
            
            for b=1:size(lbpf,1)
                xf=squeeze(lbpf_2(b,:));
                dclusfj(b)=pdist2(clusfj,xf,distanceT);
            end
            %                     xf=squeeze(lbpf(b,:));
            %                     xt=squeeze(lbpt(b,:));
            %                     if size(xf,2)~=size(clusfj,2)
            %                         if size(clusfj,2)>size(xf,2)
            %                         lbpf_N=interp1(xf,size(clusfj,2), 'spline');
            %                         else
            %                         lbpf_N=interp1(xf,size(clusfj,2), 'spline','extrap');
            %                         end
            %                         dclusfj(b)=pdist2(clusfj,lbpf_N,distanceT);
            %                     else
            %                         dclusfj(b)=pdist2(clusfj,xf,distanceT);
            %                     end
            %
            %                     if size(xt,2)~=size(clustj,2)
            %                         lbpt_N=interp1(xt,size(clustj,2), 'spline');
            %                         dclustj(b)=pdist2(clustj,lbpt_N,distanceT);
            %                     else
            %                         dclustj(b)=pdist2(clustj,xt,distanceT);
            %                     end
            %                 end
            
            dclust(j)=min(dclustj);
            dclusf(j)=min(dclusfj);
            
        end
        dref=prctile(dclust,20);
        dreff=prctile(dclusf,20);
        D(t-2,t2-2)=dref;
        Df(t-2,t2-2)=dreff;
    end
end
%end
h=figure
colormap('bone')
imagesc(max(D(:))-D)
axis off
screen_size = get(0, 'ScreenSize');
set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
if savePlots
    saveas(h,[StatsPath dataset tagfolder filesep 'past_scoreN' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
end

h2=figure
colormap('pink')
imagesc(max(Df(:))-Df)
axis off
screen_size = get(0, 'ScreenSize');
set(h2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
if savePlots
    saveas(h2,[StatsPath dataset tagfolder filesep 'future_scoreN' tagDesc '-' selCtag tag_modes '-modes' num2str(mxCluster) '.png'],'png')
end


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
D=zeros(size(clus,2), 5);
Df=zeros(size(clus,2), 5);
lag=5

%WE COMPARE AT LEAST 3 TIME STEPS
%We skip t=1 as the lagrangian is all 1

for j=1:size(clus,1)
    clustj=clus(j,:);
    
    for b=1:size(lbpRef,1)
        
        x=xcorr(clustj,lbpRef(b,:),lag,'coeff');
        %c=xcov(clustj,lbpRef(b,:),0,'coeff')
        %cc=pdist2(clustj,lbpRef(b,:),'correlation')   
        
        %xx=corrcoef(clustj,lbpRef(b,:))
        
        %r=corr(clustj',lbpRef(b,:)')
        
        %max(c)
        %max(x)
        %pause
        
        dclustj(b)=1-max(x);
        
    end
    dclust(j)=min(dclustj);
    
end

dref=prctile(dclust,25);
dref2=prctile(dclust,75);
drefmin=min(dclust);
drefmax=max(dclust);
drefm=mean(dclust);

s_yl=[0.02 0.1 0.1 0.1 0.1]

h=figure
boxplot(dclust,'colors',[0.4 0.4 0.4])

ylim([0 s_yl(c)])
set(gca,'XTick',[])
set(gca,'XTickLabel',[],'FontSize',18)
set(gca,'YTick',[0 s_yl(c)])
set(gca,'YTickLabel',[0 s_yl(c)],'FontSize',18)
            
saveas(h,[StatsPath dataset tagfolder filesep 'boxplot_score' '-' tagDesc '_' num2str(lag) '.png'],'png')
%dlmwrite([desc 'score' '-' tagDesc '_v2.csv'],D)


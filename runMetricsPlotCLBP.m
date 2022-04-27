
% /* --------------------------------------------------------------------------------------
%  * File:    runMechanomeTimeScores.m
%  * Date:    01/06/2015
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2013-2017, David Pastor Escuredo
%  with Biomedical Image Technology, UPM (BIT-UPM)
%  with BioEmergences, CNRS
%  with LifeD Lab
%  All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
dNs=[1 2 3 6 13];%first is reference


StatsPath='../MechanicsStats/';
xtag='-1-0'
useDefault=1;%Number of modes in loadStatsMetadata
derive=0
doCluster=0
doProfile=0
doBasis=1
doProjection=0
drawClusterProfile=0


tinilags=[61 62 39 8 39 95 86];%for each dataset
tinilags(13)=50;
tinilags(18)=69;
tinilags(14)=105;

%%%%%%%%%%%%%%%% Descriptor Set and Descriptor Index %%%%%%%%%%%%%%%%%%%%%
%1 Topology
%2 Strain rates
%3 Displacements Field descriptors (Unstable)
%4 Left Lagrangian descriptors
%5 Right Lagrangian descriptors
DescriptorsSets=[5];
%For 4-5
DescriptorsIndexes=[5 6 7 8];
%check libStats/loadStatsMetadata
%For 1
%DescriptorsIndexes=[4]%check libStats/loadStatsMetadata
%For 2
%DescriptorsIndexes=[1 2]%check libStats/loadStatsMetadata

close all

%%%%%%%%%%%%%%%% Tissue selection from Movit. Original domain to get modes
tagfolder='/wtmut/'
s_metric={'mean', 'prct25', 'prct75', 'min', 'max'}
s_yl=[0.003 0.007 0.004 0.005 0.005; 0.04 0.007 0.06 0.005 0.005; 0.04 0.007 0.06 0.005 0.005; 0.03 0.01 0.03 0.005 0.005]
s_yl=[0.4 0.007 0.004 0.005 0.005; 0.6 0.007 0.06 0.005 0.005; 0.6 0.007 0.06 0.005 0.005; 20 0.01 0.03 0.005 0.005]
s_mi=1

%%%%%%%%%%%%%%%%% Plotting options
extraPlots=0
seePlots=1
savePlots=1
textSimple=2
fsize=36

%%%%%%%%%%%%%%%%% time
tfinal_m=-1
tini_m=-1
hini=6
hfin=13
hpf_limits=1
crop_time_axis=1
phys_time=1
plotZero=0
plotHPF=1%Flag

step_time=20
gridOn=1
widthMode=7
use_defaultRange=1
saveClusterData=0%save
removeProfileWithInvalid=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% RUNNING %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figType=1%1 plot 2 image

globalrange=zeros(length(dNs),2)
globalrangef=zeros(length(dNs),2)

c=1
for DescriptorsSet=DescriptorsSets
    for DescriptorIndex=DescriptorsIndexes
        
        for idN=1:length(dNs)
            DescriptorIndex 
          
            dN=dNs(idN)
            tinilag=tinilags(dN)
            loadParametersKinematics
            loadStatsMetadata          
            tagDesc=['-' tagsCmaps{DescriptorIndex} tagg]
            nameDesc=tags{DescriptorIndex};
            nameDescLeg=tagsShort{DescriptorIndex};
            
            D=dlmread([desc 'CLBP_score' '-' tagDesc '_v2.csv'])
           % Df=dlmread([desc 'CLBP_future_score' '-' tagDesc '.csv'])
            
            v=D(:,s_mi);
            length(v)
            
            vmax=max(v);
            vmin=min(v);
            globalrange(idN,:)=[vmin vmax]
            h=figure
            plot(v, 'Color', [255 140 0]./255, 'LineWidth', 3)
            xlabel('time steps','FontSize',16)
            ylabel('Distance to CLBP','FontSize',16)
            ylim([0 s_yl(c,s_mi)])
            
            set(gca,'XTick',[0:20:length(v)])
            set(gca,'XTickLabel',[0:20:length(v)],'FontSize',18)
            %title([dataset tagDesc])
            grid on
            [StatsPath dataset tagfolder filesep 'plot_CLBP_score' '-' tagDesc '-' s_metric(s_mi) '.png']
            saveas(h,[StatsPath dataset tagfolder filesep 'plot_CLBP_score' '-' tagDesc '-' s_metric{s_mi} '.png'],'png')
            
% %             v2=Df(:,s_mi);
% %             length(v2)
% %             
% %             vmax=max(v2);
% %             vmin=min(v2);
% %             globalrangef(idN,:)=[vmin vmax]
% %             h=figure
% %             plot(v2, 'Color', [0 0 128]./255, 'LineWidth', 3)
% %             xlabel('time steps','FontSize',16)
% %             ylabel('cosine distance','FontSize',16)
% %             ylim([0 s_yl(c,s_mi)])
% %             
% %             set(gca,'XTick',[0:20:length(v2)])
% %             set(gca,'XTickLabel',[0:20:length(v2)],'FontSize',18)
% %             title([dataset tagDesc])
% %             grid on
% %             [StatsPath dataset tagfolder filesep 'plot_CLBP_future_score' '-' tagDesc '-' s_metric(s_mi) '.png']
% %             saveas(h,[StatsPath dataset tagfolder filesep 'plot_CLPB_future_score' '-' tagDesc '-' s_metric{s_mi} '.png'],'png')
% %             
        end  
        c=c+1
    end
end

globalrange

% maxsD=[0 0 0 0 0.0002 0.006 0.006 0.007 0]
% max(D(:))
% max(Df(:))
% 
% h=figure
% colormap('bone')
%Df=Df./maxsD(DescriptorIndex);

% if normType==1
%     
%     D=D./maxsD(DescriptorIndex);
%     D=1-D;
%     tagn='1'
% else
%     'hols'
% end

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
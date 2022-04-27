% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draws the profile of domains defined by selections %%%s
clear m
clear m3
%%%% Make the computations if necessary
if ~exist([stats_file '.mat'])
   'NO EXIST'
   %pause
   doProfile=1
   loadStatsMetadata
   runDescriptorsProfile 
end
%m_title=[m_title '-' dd]
loadPlotConfiguration;

%%%% Load metric data structure
load([stats_file '.mat'])

%%%% Time scale
%discrete_time
%distime=[0:size(metric,1)-1];
distime=[tini:tini+size(metric,1)-1]
if T>1 && downsampleT
    sampl=[1:T:size(metric,1)]
    metric=metric(sampl,:,:);
    t_start=floor(t_start/T)+1
    t_limit=floor(t_limit/T)+1   
end

%%%% Setting timeline
phystime=(distime*dt_min)+t0_min;
if phys_time
    time_axis=phystime;
else                         
    time_axis=distime;
end
max(time_axis)
min(time_axis)

%%%% Select data
%Time %Stat % Selection
size(metric)
%the statSelector should have the main stat in 1
m=metric(:,statSelector,:);

%average 2 populations
if length(aveSelector)==2
    m(:,:,aveSelector(1))=(m(:,:,aveSelector(1))+m(:,:,aveSelector(2)))/2;
    m(:,:,aveSelector(2))=[];
end
m=m(:,:,selSelector);
size(m)

%%%% Crop plots in time (still the original axis)
if hpf_limits
   loadTimelineConfiguration 
end
tini_mat=t_start+1
tfin_mat=t_limit+1
tinilag_mat=tinilag+1
%zv=m(tini_mat,1,1);%maybe not right
m3=zeros(size(m));%zero padding 

m3(tini_mat:tfin_mat,:,:)=m(tini_mat:tfin_mat,:,:);
samples=metric(tini_mat:tfin_mat,12,selSelector);
if DescriptorsSet>3 
    m3(1:tinilag_mat-1,:,:)=NaN;
    m3(tinilag_mat,1,:)=zv(DescriptorIndex);
end
    
%Boundaries of the data in y axis
data=m3(tini_mat:tfin_mat,plotStatSelector,:);
maxdata=nanmax(data(:))%limits of the data
mindata=nanmin(data(:))%limits of the data

%Not just the plot but the whole figur is cropped in time
if crop_time_axis
    m3=m3(tini_mat:tfin_mat,:,:);
    if phys_time
        time_axis=time_axis(tini_mat:tfin_mat);
    else
        time_axis=time_axis(tini_mat:tfin_mat);
    end
end

if startInTref
    time_axis=[1:1:length(time_axis)]
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING DESCRIPTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%in case we want to overlay different plots
if ~overlayed
    h=figure
    hold on
end
if hidePlot
    set(h,'Visible','off')
end

for ss=1:length(plotSelector)
    for sst=1:length(plotStatSelector)
        sp=plotSelector(ss);
        %%%%%% restrict min samples
        samplessel=squeeze(samples(:,1,sp));
        m3(samplessel<min_samples,:,sp)=NaN;
        %%%%% color
        m_color=cl(sp,:)+(sst-1)*[0.35 0.3 -0.2];
        m_color(m_color>1)=1;
        m_color(m_color<0)=0;
        %plot(m3(:,ss),'Color',[1 0 0],'LineWidth',2)
        if length(errorSelector)>0
            shadedErrorBar(time_axis, squeeze(m3(:,plotStatSelector(sst),sp)),squeeze(m3(:,errorSelector(sst),sp)),{'Color',m_color,'LineWidth',tl},1)%,),'LineWidth',tl)
        else
            plot(time_axis, squeeze(m3(:,plotStatSelector(sst),sp)),'Color',m_color,'LineWidth',tl)
        end
    end
end

%Legend
if simple<2
legend(my_legend)
end

marginM=0.01
marginm=-0.01;
if DescriptorsSet>3 
    marginM=0.1
    marginm=-0.1;  
end

mindata=desc_limits(DescriptorIndex,1)
maxdata=desc_limits(DescriptorIndex,2)

if mindata==-1 && maxdata==-1
    desc_limits(DescriptorIndex,:)=[mindata+marginm maxdata+marginM]
end

%%%% X (TIME AXIS) CONFIGURATION!!!
if phys_time    
    if simple<2
    xlabel('Time (mins)','FontSize', fsize)
    set(gca,'XTick',time_axis([1:10:end]))
    set(gca,'XTickLabel',time_axis([1:10:end]),'FontSize',fsize) 
    end
    %set(gca,'XTickLabel',{'hpf7','hpf8','hpf9','hpf10','hpf11','hpf12'}) 
    if plotHPF
        plot([hpf7 hpf7],[mindata maxdata], '--','Color',hc,'LineWidth',hl);
        plot([hpf8 hpf8],[mindata maxdata], '--','Color',hc,'LineWidth',hl);
        plot([hpf9 hpf9],[mindata maxdata],'--','Color',hc,'LineWidth',hl);
        plot([hpf10 hpf10],[mindata maxdata], '--','Color',hc,'LineWidth',hl);
        plot([hpf11 hpf11],[mindata maxdata], '--','Color',hc,'LineWidth',hl);
        plot([hpf12 hpf12],[mindata maxdata], '--','Color',hc,'LineWidth',hl);
        plot([hpf13 hpf13],[mindata maxdata], '--','Color',hc,'LineWidth',hl);
        plot([hpf14 hpf14],[mindata maxdata], '--','Color',hc,'LineWidth',hl);
        set(gca,'XTick',[hpf6 hpf7 hpf8 hpf9 hpf10 hpf11 hpf12 hpf13 hpf14])
        set(gca,'XTickLabel',{'6', '7','8','9','10','11','12','13','14'},'FontSize',fsize) 
        set(gca,'XLim',[time_axis(1) hpf14+1]) 
        if simple<2
            xlabel('Time (hpf)','FontSize',fsize)
        end
    end
    
    if printTrefs
       if startInTref
        plot([1 1].*dt_min+t0_min, [mindata maxdata], '--','Color',[0 0 0],'LineWidth',hl+1); 
        set(gca,'XTick',[1])
        set(gca,'XTickLabel',{'tini'},'FontSize',fsize)
       else
       plot([tinilags(dN) tinilags(dN)].*dt_min+t0_min, [mindata maxdata], '--','Color',[0 0 0],'LineWidth',hl+1);
       end
    end
else
    if printTrefs
        if startInTref
        plot([1 1], [mindata maxdata], '--','Color',[0 0 0],'LineWidth',hl+1); 
        set(gca,'XTick',[])
       set(gca,'XTickLabel',{' '},'FontSize',fsize)
       % set(gca, 'TickLength', [0 0]);
        for lx=[1:5]
        taa=floor(lx*60/2.5)
        plot([taa taa], [mindata maxdata], '--','Color',hc,'LineWidth',hl); 
       % set(gca,'XTick',[1 taa])
        end
        %set(gca,'XTickLabel',{'tini' 'tini+5'},'FontSize',fsize)
        %set(gca, 'interpreter', 'tex');
        else
            plot([tinilags(dN) tinilags(dN)], [mindata maxdata], '--','Color',[0 0 0],'LineWidth',hl+1);
        end
    end
    if simple<2
    xlabel('Time Step','FontSize',fsize) 
    set(gca,'XTick',time_axis([1:10:end]))
    set(gca,'XTickLabel',time_axis([1:10:end])) 
    end
end
if simple==3 && ~startInTref
    set(gca,'xtick',[],'FontSize',fsize)
end

%%%% Y AXIS CONFIGURATION! 
%ylabel(m_ylabel,'FontSize',15)
ylim(desc_limits(DescriptorIndex,:))
set(gca,'YTick',desc_limits(DescriptorIndex,:))
ymnlabel=desc_limits(DescriptorIndex,1)
ymxlabel=desc_limits(DescriptorIndex,2)
if ymnlabel<0.01
    ymnlabel = sprintf('%.g',desc_limits(DescriptorIndex,1));
else
    ymnlabel = num2str(desc_limits(DescriptorIndex,1));
end
if ymxlabel<0.01
    ymxlabel = sprintf('%.g',desc_limits(DescriptorIndex,2));
else
    ymxlabel = num2str(desc_limits(DescriptorIndex,2));
end

if setYaxis
    if writeZero && (desc_limits(DescriptorIndex,1) <zv(DescriptorIndex))     
        set(gca,'YTick',[desc_limits(DescriptorIndex,1) zv(DescriptorIndex) desc_limits(DescriptorIndex,2)]) 
        set(gca,'YTickLabel',{ymnlabel, num2str(zv(DescriptorIndex)), ymxlabel}, 'FontSize',fsize)
    else
        set(gca,'YTick',desc_limits(DescriptorIndex,:) ) 
        set(gca,'YTickLabel',{ymnlabel, ymxlabel}, 'FontSize',fsize) 
    end
end

if plotZero
    zl=0.5
    if phys_time
        plot([hini*60 hfin*60],[zv(DescriptorIndex) zv(DescriptorIndex)], '--','Color','k','LineWidth',zl);
    else
        plot([0 t_limit-t_start],[zv(DescriptorIndex) zv(DescriptorIndex)], '--','Color','k','LineWidth',zl);
    end
end

%%%% TEXT AND GENERAL
if simple<2
    title([dataset ' ' m_title],'FontSize',fsize+4)
    m_title
    if simple==1
       title(['Dataset ' num2str(dN) ' - ' tagDescOr],'FontSize',fsize+4) 
    end
end

if gridOn
    grid On
    set(gca,'GridLineStyle',gridStyle)
else
    grid Off
end

if ~overlayed
    hold off
end
%%%% Screen size
screen_size = get(0, 'ScreenSize');
set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%%%% Paper size
            %r=120%resolution screen-paper
            %set(h, 'PaperUnits','inches','PaperPosition', [0 0 screen_size(3) screen_size(4) ]/r );
%%%% Save plot
if savePlot && ~overlayed
    saveas(h,[profile_file '.png'],'png');
    %print(h,'-dpng',[m_title '.png'])
end


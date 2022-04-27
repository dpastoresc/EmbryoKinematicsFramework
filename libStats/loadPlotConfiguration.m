% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Plots configuration %%%%%%%%%%%%%%%%%%%%%

%%%% Defaults
m_title='Descriptor kinematics'
m_ylabel='Magnitude axis'
m_legend=''

%%%% Configurations
% Colors
cl=['m','b','k','g'];
cl=['b','g','r','m'];
cl=[1 0 0 ; 0 0 1;  0 1 0; 0.5 0.6 0.7]
cl=[255 218 88; 0.1 1 0.1; 0.2 0.5 0.98; 0.5 0.6 0.7]
cl=[200 190 88; 0.1 1 0.1; 0.2 0.5 0.98; 0.5 0.6 0.7]
cl(1,:)=cl(1,:)/255.0
cl=['r','b','r','m'];
cl=[1 0 0 ; 0 1 0;  0 0 1; 0.5 0.6 0.7]
cl=[200 154 0 ; 0 0.7 0;  0 0 1; 0.5 0.6 0.7]
cl(1,:)=cl(1,:)/255.0
%cl=[0.7 0.2 0.8; 0 1 0; 0.5 0.3 1; 0 0 1]
if drawClusterProfile
    cl=[0 0 0.8; 0 1 0; 0.4 0.3 0; 1 0 1; 0.6 0 0]
end
tl=2.5

%%%% Selector of data. Default ones
aveSelector=[]%stats to be averaged
statSelector=[1]%stats u keep
errorSelector=[]%stats used to plot the error profile (1 symmetric 2 assymetric)
plotStatSelector=[1]%stats to be ploted wihting statSelector
%initialValueFactor=1;

tissues2={'Hypoblast', 'Eye', 'Epiblast', 'All', 'Midline'}

%%%% Stats
%stat={'Mean','Median','std','Max',',Min','Percentile-99','Percentile-1' ...
%    'Percentile-75','Percentile-25','entries','entries-dif','samplesTaken' ...
%    'samplesInside', 'ss', 't'} 

%%%% Configurations setups
%Plotting simple mean median
if config==1
    %selections
    selSelector=[1 3]
    plotSelector=[1 2]%relative to selSelector
    m_legend=tissues2; 
    %configure stast
    statSelector=[1 2]%relative to statSelector
    errorSelector=[]%relative to statSelector
    plotStatSelector=[1 2]%relative to statSelector
    %names and colors
    cl=cl([1 3],:)
    m_title=stats_options%[dd ' ' tagDesc]

elseif config==11
    %selections
    selSelector=[1 3]
    plotSelector=[1 2]%relative to selSelector
    m_legend=tissues2; 
    %configure stast
    statSelector=[1 3]%relative to statSelector
    errorSelector=[]%relative to statSelector
    plotStatSelector=[1 2]%relative to statSelector
    %names and colors
    cl=cl([1 3],:)
    m_title=stats_options%[dd ' ' tagDesc]
    desc_limits=[-1 -1]
    
elseif config==13
    %selections
    selSelector=[1 3]
    plotSelector=[1 2]%relative to selSelector
    m_legend=tissues2; 
    %configure stast
    statSelector=[3]%relative to statSelector
    errorSelector=[]%relative to statSelector
    plotStatSelector=[1]%relative to statSelector
    %names and colors
    cl=cl([1 3],:)
    m_title=stats_options%[dd ' ' tagDesc]
    desc_limits=[-1 -1]
    
    
%Plotting with errors std
elseif config==2
    %selections
    selSelector=[1 3]
    plotSelector=[1 2]%relative to selSelector
    %configure stast
    statSelector=[1 3]%relative to statSelector
    errorSelector=[2]%relative to statSelector
    %names and colors
    m_legend=tissues2; 
    cl=cl([1 3],:)
    m_title=stats_options%[dd ' ' tagDesc]  
    
    %Plotting with errors std
    %Plotting with errors std
elseif config==222
    %selections
    selSelector=[2 3]
    plotSelector=[1 2]%relative to selSelector
    %configure stast
    statSelector=[1 3]%relative to statSelector
    errorSelector=[2]%relative to statSelector
    %names and colors
    m_legend=tissues2; 
    cl=cl([2 3],:)
    m_title=stats_options%[dd ' ' tagDesc]  
    
    %Plotting with errors std
elseif config==21
    %selections
    selSelector=[3]
    plotSelector=[1]%relative to selSelector
    %configure stast
    statSelector=[1 3]%relative to statSelector
    errorSelector=[2]%relative to statSelector
    
    %names and colors
    m_legend=tissues1; 
    cl=cl([3],:)
    m_title=stats_options%[dd ' ' tagDesc] 
    
        %Plotting with errors std
elseif config==23
    %selections
    selSelector=[3]
    plotSelector=[1]%relative to selSelector
    %configure stast
    statSelector=[3]%relative to statSelector
    errorSelector=[]%relative to statSelector
    plotStatSelector=[1]%relative to statSelector
    %names and colors
    m_legend=tissues1; 
    cl=cl([3],:)
    m_title=stats_options%[dd ' ' tagDesc] 
    
%Plotting with errors 25-75 percentiles
elseif config==3
    %selections
    selSelector=[1 3]
    plotSelector=[1 2]%relative to selSelector
    %configure stast
    statSelector=[1 8 9]
    errorSelector=[2 3]%relative to statSelector
    %names and colors
    m_legend=tissues2; 
    cl=cl([1 3],:)
    m_title=stats_options%[dd ' ' tagDesc]  
    
%Plotting samples
elseif config==4
    %selections
    selSelector=[1]
    plotSelector=[1]%relative to selSelector
    %configure stast   
    statSelector=[10:14]%relative to statSelector
    plotStatSelector=[ 1 3 4]%relative to statSelector
    %names and colors
    cl=cl([1 3],:)
    m_legend=tissues2; 
    m_title=stats_options%[dd ' ' tagDesc]
    desc_limits=-1.*ones(8,2);
    
elseif config==99
    %selections
    selSelector=[1:mxCluster]
    plotSelector=[1:mxCluster]%relative to selSelector
    %configure stast   

    statSelector=[1 3]%relative to statSelector
    errorSelector=[2]%relative to statSelector
    %names and colors
    %names and colors
    cl=cl([1:mxCluster],:)
    for ixx=1:mxCluster
        m_legend{ixx}=['Class ' num2str(ixx)]
    end        
    m_title=stats_options%[dd ' ' tagDesc]
    
elseif config==50
    %selections
    selSelector=[10:12]
    plotSelector=[1:3]%relative to selSelector
    %configure stast   

    statSelector=[10:14]%relative to statSelector
    plotStatSelector=[2 3 4]%relative to statSelector
    %names and colors
    cl=cl([1 3],:)
    m_legend=tissues2; 
    m_title=stats_options%[dd ' ' tagDesc]
end

%%%% Legend names configuration
cleg=1
if length(errorSelector)>0
    for l=1:length(m_legend)
        my_legend{cleg}=[m_legend{l}]              
        if length(errorSelector)>1
            my_legend{cleg+1}=stat{statSelector(errorSelector(1))}
            my_legend{cleg+2}=stat{statSelector(errorSelector(2))}
        else
            my_legend{cleg+1}=['+ ' stat{statSelector(errorSelector(1))}]
            my_legend{cleg+2}=['- ' stat{statSelector(errorSelector(1))}]
        end
        my_legend{cleg+3}=stat{statSelector(1)}  
        cleg=cleg+4
    end
else
    for l=1:length(m_legend)
        for sll=1:length(plotStatSelector)
            my_legend{cleg}=[m_legend{l} ' ' stat{statSelector(plotStatSelector(sll))}]
            cleg=cleg+1
        end
    end
end



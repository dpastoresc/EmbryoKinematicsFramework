% /* --------------------------------------------------------------------------------------
%  * File:    datasetListBioEmergences.m
%  * Date:    01/06/2015
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2013-2017, David Pastor Escuredo
%  with Biomedical Image Technology, UPM (BIT-UPM)
%  with BioEmergences, CNRS
%  All rights reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ZEBRAFISH BIOEMERGENCES DATASETS INPUT
%Set your dN in loadDataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataLocal=0;
dataNewOrientation=0
newData=0

if dataNewOrientation
    if dN==2
        dataset='151113a'
        trackID='21048'
        tini=4
        tfinal=131
        if tinilag<0
            tinilag=0
        end
    elseif dN==1
        dataset='151113aM'
        trackID='21134'
        tini=4
        tfinal=108
        if tinilag<0
            tinilag=0
        end
    elseif dN==3
        dataset='151118aM'
        trackID='21135'
        tini=10
        tfinal=139
        if tinilag<0
            tinilag=0
        end
    elseif dN==4
        dataset='151124aF'%not found
        trackID='21137'
        tini=28
        tfinal=150
        if tinilag<0
            tinilag=0
        end
    elseif dN==5
        dataset='151125aM'
        trackID='21137'
        tini=4
        tfinal=120
        if tinilag<0
            tinilag=0
        end
    elseif dN==6
        dataset='151126aF'%not found
        trackID='21138'
        tini=52
        tfinal=173
        if tinilag<0
            tinilag=0
        end
    elseif dN==7
        dataset='160524aF'
        trackID='21168'
        tini=52
        tfinal=173
        if tinilag<0
            tinilag=0
        end
    elseif dN==8
        dataset='160525aF'
        trackID='21167'
        tini=51
        tfinal=171
        if tinilag<0
            tinilag=0
        end
        
    elseif dN==9
        dataset='160105aF'
        trackID='21133'
        tini=34
        tfinal=150
        if tinilag<0
            tinilag=0
        end
        phenotype='trilobite'
        
    elseif dN==10
        dataset='160105a'
        trackID='21139'
        tini=37
        tfinal=150
        if tinilag<0
            tinilag=0
        end
        phenotype='trilobite'
        
    elseif dN==11
        dataset='170123a'
        trackID='21233'
        tini=0
        tfinal=350
        if tinilag<0
            tinilag=0
        end
        phenotype='syamese'
    end
else

    if dN==1
        dataset='071226a'
        trackID='20801'
        tini=0;
        tfinal=230;
        if tinilag<0
            tinilag=108
            %61 71 85 108 118 131 155
            % 178 189 201
        end
        t_max_cluster=178
%         t_max_cluster=85
%         t_max_cluster=108
%         t_max_cluster=131
%         t_max_cluster=155
%         t_max_cluster=201
    elseif dN==111%tracking alternative 1
        dataset='071226a'
        trackID = '20046'
        tini=0;
        tfinal=230;
        if tinilag<0
            tinilag=51
        end
        
    elseif dN==112%tracking alternative 2
        dataset='071226a'
        trackID = '20836'
        tini=0;
        tfinal=230;
        if tinilag<0
            tinilag=51
        end
        
    elseif dN==2
        dataset='140523aF'
        trackID='20787'
        tini=0;
        tfinal=209;
        if tinilag<0
            tinilag=62
        end
        t_max_cluster=209
        
    elseif dN==3
        dataset='141106aF'
        trackID='20849'
        tini=0
        tfinal=200%230
        if tinilag<0
            tinilag=39
        end
        t_max_cluster=200
        
    elseif dN==4
        dataset='141108aF'
        trackID='20854'
        tini=0
        tfinal=168%230
        if tinilag<0
            tinilag=8
        end
        t_max_cluster=168
        
    elseif dN==5
        dataset='141108a'
        trackID='20865'
        tini=0
        tfinal=198%230
        if tinilag<0
            tinilag=39
        end
        t_max_cluster=198
        
    elseif dN==91
        dataset='141121a'
        trackID='20886'
        tini=0
        tfinal=265
        if tinilag<0
            tinilag=0
        end
        
    elseif dN==92
        dataset='141121aF'
        trackID='20888'
        tini=0
        tfinal=265
        if tinilag<0
            tinilag=0
        end
        
    elseif dN==63
        dataset='081018a'
        trackID='21188'
        tini=0
        tfinal=237%200
        if tinilag<0
            tinilag=106
        end
        dataLocal=0;
        
    elseif dN==61
        dataset='081018a'
        trackID='20055'
        tini=0
        tfinal=237%200
        if tinilag<0
            tinilag=106
        end
        dataLocal=0;
        
    elseif dN==62
        dataset='081018a'
        trackID='20972'
        tini=0
        tfinal=237%200
        if tinilag<0
            tinilag=106
        end
        dataLocal=0;
        
    elseif dN==6
        dataset='081018a'
        trackID='21189'
        tini=0
        tfinal=237%200
        if tinilag<0
            tinilag=154
            %95 106 130 154 178
        end
        t_max_cluster=213%13hpf
        dataLocal=0;
        
    elseif dN==71
        dataset='081025a'
        trackID='20056'
        tini=0
        tfinal=320%200
        if tinilag<0
            tinilag=86
        end
        dataLocal=0;
        
    elseif dN==7
        dataset='081025a'
        trackID='21211'
        tini=0
        tfinal=222%320%
        if tinilag<0
            tinilag=154
            %86, 108, 154
        end
        dataLocal=0;
        t_max_cluster=222%13hpf
        
    elseif dN==72
        dataset='081025a'
        trackID='21221'
        tini=0
        tfinal=320%200
        if tinilag<0
            tinilag=86
        end
        dataLocal=0;
        
    elseif dN==8
        dataset='160330aF'
        trackID='21199'
        tini=0
        tfinal=222%200
        if tinilag<0
            tinilag=50
        end
        dataLocal=0;
        
    elseif dN==81
        dataset='160330aF'
        trackID='21201'
        tini=0
        tfinal=222%200
        if tinilag<0
            tinilag=50
        end
        dataLocal=0;
        
    elseif dN==10
        dataset='160324aF'
        trackID='21191'
        tini=0
        tfinal=350%200
        if tinilag<0
            tinilag=48
        end
        dataLocal=0;
        
    elseif dN==101
        dataset='160324aF'
        trackID='21192'
        tini=0
        tfinal=350%200
        if tinilag<0
            tinilag=48
        end
        dataLocal=0;
        
    elseif dN==11
        dataset='161012a'
        trackID='21197'
        tini=0
        tfinal=255%200
        if tinilag<0
            tinilag=57
        end
        dataLocal=0;
        
    elseif dN==12
        dataset='100210aF'
        trackID='20295'
        tini=0
        tfinal=225%200
        if tinilag<0
            tinilag=45
        end
        dataLocal=0;
        
    elseif dN==121
        dataset='100210aF'
        trackID='20458'
        tini=0
        tfinal=225%200
        if tinilag<0
            tinilag=45
        end
        dataLocal=0;
        
    elseif dN==13
        dataset='120914aF'
        trackID='20484'
        tini=0
        tfinal=169%193%200
        if tinilag<0
            tinilag=109
            %8-50  %10.5 109
        end
        dataLocal=0;
        t_max_cluster=169 %13hpf
        
    elseif dN==14
        dataset='091021aF'
        trackID='20725'
        tini=0
        tfinal=184%204
        if tinilag<0
            tinilag=105
            %105 and 154 (10.5) and 184 (12)
        end
        dataLocal=0;
        t_max_cluster=200%there is something weird at the end so we save 4 steps
        t_max_cluster=184
        tfinal=184
    elseif dN==15
        dataset='091016aF'
        trackID='21229'
        tini=0
        tfinal=213%200
        if tinilag<0
            tinilag=154%103%154
            %103 8 %158 10.5
        end
        dataLocal=0;
        t_max_cluster=208
        
    elseif dN==151
        dataset='091016aF'
        trackID='20121'
        tini=0
        tfinal=270%200
        if tinilag<0
            tinilag=45
        end
        dataLocal=0;
        
    elseif dN==16
        dataset='170315aF'
        trackID='21252'
        tini=0
        tfinal=250%286%250(13)
        if tinilag<0
            tinilag=73%73%143
        end
        dataLocal=0;
        
    elseif dN==17
        dataset='170321aF'
        trackID='21250'
        tini=0
        tfinal=270%243, 216
        if tinilag<0
            tinilag=135%81%135
        end
        dataLocal=0;
    
    elseif dN==18
        dataset='071224aF'
        trackID='21272'
        tini=0
        tfinal=172%155
        if tinilag<0
            tinilag=112
            %69 %112
        end
        dataLocal=0;
        t_max_cluster=155
   
    elseif dN==19
        dataset='091023aF'
        trackID='21274'
        tini=0
        tfinal=172%243, 216
        if tinilag<0
            tinilag=69%
        end
        dataLocal=0;
        
    elseif dN==20
        dataset='160106'
        trackID=''
        tini=0
        tfinal=49%243, 216
        if tinilag<0
            tinilag=0%81%135
        end
        dataLocal=1;
    end
    
end

if newData
    if dN==1
        dataset='180516aF'
        trackID='21383'
        tini=0
        tfinal=300
        if tinilag<0
            tinilag=0
        end
    elseif dN==2
        dataset='180420hZ'
        trackID='21361'
        tini=0
        tfinal=300
        if tinilag<0
            tinilag=0
        end
    end
end








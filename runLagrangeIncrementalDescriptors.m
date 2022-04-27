% /* --------------------------------------------------------------------------------------
%  * File:    runLangrangeIncrementalDescriptors2.m
%  * Date:    01/06/2015
%  * Author:  David Pastor Escuredo, david.pastor@upm.es
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015, David Pastor Escuredo - david.pastor@upm.es
%  with Biomedical Image Technology, UPM (BIT-UPM)
%  All rights reserved.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAGRANGE KINEMATICS DESCRIPTORS FROM TRACKING
% This version incrementes the temporal window of acummulation T(t) from a reference point

clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dNs=[14] %SET DATASET ID (datasetListBioEmergences.m)
process=1%just process (1) or export to movit (0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for dN=dNs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tinilag=-1%SET or use default (-1) in datasetListBioEmergences.m
    real_pos=1
    reuse_links=0
    xtag=['-' num2str(real_pos) '-' num2str(reuse_links)]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    loadParametersKinematics;
    tinilag
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% READ ALL THE INTERMEDIATE DATA %%%%%%%%%%%%
    'Read Deformation Tensor Data'
    
    F=load([folder dataset '_' tensorid '.mat']);
    F = struct2cell(F);
    F=F{1};
    XYZ=load([folder dataset '_XYZ.mat']);
    XYZ = struct2cell(XYZ);
    XYZ=XYZ{1};
    cellid=load([folder dataset '_cellid.mat']);
    cellid = struct2cell(cellid);
    cellid=cellid{1};
    %cellnum=load_centData(dataset,'cellnum');
    cellnum=load([folder dataset '_cellnum.mat']);
    cellnum = struct2cell(cellnum);
    cellnum=cellnum{1};
    cellmother=load([folder dataset '_mother.mat']);
    cellmother = struct2cell(cellmother);
    cellmother=cellmother{1};
    cellchild=load([folder dataset '_child.mat']);
    cellchild = struct2cell(cellchild);
    cellchild=cellchild{1};
    Velo = load([folder dataset '_' veloid]);%load_centData(dataset, veloid);
    Velo = struct2cell(Velo);
    Velo=Velo{1};
    MeanVelo=load([folder dataset '_' velomeanid]);%load_centData(dataset, velomeanid);
    MeanVelo = struct2cell(MeanVelo);
    MeanVelo=MeanVelo{1};
    N=load([folder dataset '_' neighid]);%load_centData(dataset, neighid);
    N = struct2cell(N);
    N=N{1};
    
    %Initialization of data structure
    nl=33;
    stepcheck=32;%number of Ts accumulated so the integration is coherent
    
    %identity matrix
    I=[1 0 0; 0 1 0; 0 0 1];
    dt_step=0;
    if update_pos
        dt_step=dt_min;
    end
    DescriptorsL=[];
    DescriptorsLLeft=[];
    
    %Control Variables
    countNullTensorTotal=0
    debug_this=0
    
    %matlab indexes are 1-based
    tini_mat=tinilag+1;
    
    %Save the dynamic lagrngian tracking from tinilag
    if save_tracking
        % Tracking format (.csv)
        % track_id;cellid;time_step;
        track_file=[desc 'tracking-' num2str(tinilag) '-' num2str(tfinal) '-' num2str(real_pos) '-' num2str(reuse_links) '.csv']
        %cellid;trackid;t
        track_lag=double([])
    end
    
    if process
        %FIRST LOOP!
        %We update tinilag->tini2 T(t)
        for tcurrent=tinilag:step:tfinal
            'Accumulation from tinilag to current step'
            tcurrent
            %tcurrent=tinilag->Tt=0 no accumulation -> Identity
            Tt=tcurrent-tinilag
            tcurrent_mat=tcurrent+1;%matlab 1-based
            
            %from t1->t0, we use the tensors of both steps
            tlimit_track=tinilag+Tt
            
            %Read origin and destination coordinates
            coords=XYZ{tini_mat};
            cells=size(coords,2);
            coordsEnd=XYZ{tcurrent_mat};
            cellsEnd=size(coordsEnd,2);
            
            %Init data structures
            DescriptorsLs=double(zeros(cells-1,nl));
            DescriptorsLsLeft=double(zeros(cells-1,nl));
            count=1;
            
            %MATERIAL ID / CELLID LOOP!
            %We create a pathline for every cell position in tinilag
            %First cell is nothing, control row, so we iterate from pos 2
            for cell=2:cells
                %Each cell in the first frame starts a track_id
                track_id=cell-1;%starting in 1
                
                %%%%%%%%%%%%%%% The tensor in the original timestep is identity %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %init kin datas
                globalTensor=I;
                lagvelo=0;
                acP=0;
                countGoodSteps=0;%control counts
                countNullTensor=0;
                
                %%%%%%%%%%%%%%%%%%%%% Initialize building the propagation %%%%%%%%%
                %Initiates the cell origin index in tinilag of the current
                %tracking arrays
                cell_ini=cell;
                %Position of the mother cell to have the deformation to the child
                cellprev=cell;
                
                if save_tracking && tcurrent==tfinal
                    track_cells_ids=zeros(Tt,8);
                    track_cells_ids(1,:)=[track_id cellid{tini_mat}(:,cell_ini) XYZ{tini_mat}(:,cell_ini)' tinilag 0 cell_ini];
                    %track_cells_ids2=zeros(Tt,8);
                    %track_cells_ids2(1,:)=[track_id cellid{tini_mat}(:,cell_ini) XYZ{tini_mat}(:,cell_ini) tinilag 0 cell_ini];
                    count_track=2;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUILD COMPLETE TRAJECTORY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %INTERNAL TRAJECTORY LOOP!!
                %For T=0-> t_destiny=inilag-> no tracking->no accumulation->identity deformation
                
                %The position in t uses the accumulation from tinilag to t-1 as in
                %t-1 we have the deformation from t-1 -> t
                for t_destiny=(tinilag+1):tcurrent
                    
                    %'accumulating'
                    t_des_mat=t_destiny+1;%matlab is 1-based
                    t_prev_mat=t_des_mat-1;%the deformation from t-1->t is expressed in t
                    
                    %%%%%%%%%%%%%%%% Get the deformation tensor in t
                    locTensor=zeros([3 3]);
                    %Def data in t-1->t in the right index
                    locTensor(:)=F{t_prev_mat}(:,cellprev);
                    
                    if ~any(locTensor~=0)
                        %'ALL ZEROS'
                        %we lose a lot of cells when G is small
                        countNullTensorTotal=countNullTensorTotal+1;
                        countNullTensor=countNullTensor+1;
                        locTensor=I;
                    else
                        %Chain rule
                        globalTensor=locTensor*globalTensor;
                        %count
                        countGoodSteps=countGoodSteps+1;
                        %accumulated P. Just to check
                        acP=acP+abs(trace(locTensor-I));
                    end
                    
                    %%%%%%%%%%%%%%% Find the current position %%%%%%%%%%%%%%%%%%%%%
                    %Get data for the current cell in this loop
                    childIndex=cellchild{t_prev_mat}(:,cellprev);
                    childIndex=childIndex(childIndex>1);
                    if real_pos
                        %Closest neighbour interpolation
                        if t_destiny>(tinilag+1)
                            my_pos=next_pos;
                            if ~reuse_links
                                %we can keep the cell tracking links or rely in the vector
                                %field integration so all links are type 2
                                childIndex=[];
                            end
                        else
                            my_pos=XYZ{t_prev_mat}(:,cellprev);
                        end
                    else
                        %Only samples within the cell tracking
                        my_pos=XYZ{t_prev_mat}(:,cellprev);
                    end
                    
                    %The trajectory interpolates the information using the
                    %regularized flow vector field
                    my_vel=MeanVelo{t_prev_mat}(1:3,cellprev);
                    displace=my_vel*dt_min;
                    next_pos=my_pos+displace;
                    %data
                    lagvelo=norm(my_vel)+lagvelo;
                    link=0;
                    
                    %%%%%%%%%%%% Have the closest cell in the tracking %%%%%%%%%%%%
                    %%%%%%%%%%%%% Closest tracked cell interpolation %%%%%%%%%%%%%%
                    if size(childIndex,1)>1                        
                        celXYZ=XYZ{t_des_mat}(:,childIndex);
                        celvel=MeanVelo{t_des_mat}(1:3,childIndex);
                        %best child according to distance and vectory similarity
                        ca=[acosd(dot(celvel(:,1)./norm(celvel(:,1)),my_vel./norm(my_vel))),acosd(dot(celvel(:,2)./norm(celvel(:,2)),my_vel./norm(my_vel)))];
                        da=[norm(celXYZ(:,1)-(my_pos+displace)), norm(celXYZ(:,2)-(my_pos+displace))];
                        if abs(diff(da))>5 || abs(diff(ca))<3
                            bestChild=find(da==min(da));
                        else
                            bestChild=find(ca==min(ca));
                        end
                        childIndex=childIndex(bestChild(1));
                        link=1;
                    end
                    
                    %This means no children is found. Thee tracking original is cut
                    %so we need to create a link dinamically
                    if size(childIndex,1)==0
                        next_cand=XYZ{t_des_mat};
                        next_cand_vel=MeanVelo{t_des_mat};
                        %COMPLETING TRAJECTORY: mechanism to chose the next pos in the trajectory
                        %dinamically
                        if weight_similar_vel
                            %We consider also the similarity of velocity and not
                            %just the distance
                            childIndex=get_dynamic_child(my_pos, my_vel, next_cand, next_cand_vel,dt_step);
                        else
                            childIndex=get_dynamic_child(my_pos, next_cand, my_vel, dt_step);%updating position
                        end
                        link=2;
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if save_tracking && tcurrent==tfinal
                        %insert the cell id in the tracking data
                        m_cellid=cellid{t_des_mat}(:,childIndex);
                        track_cells_ids(count_track,[1 2 6 7 8])= double([track_id m_cellid t_destiny link childIndex]);
                        track_cells_ids(count_track,3:5)=double(next_pos');
                        count_track=count_track+1;
                    end
                    cellprev=childIndex;    
                end
                
                %put all trajectories together one after the other
                if save_tracking && tcurrent==tfinal
                    track_lag=vertcat(track_lag, track_cells_ids);
                    track_id=track_id+1;
                end
                
                [Lag, LagLeft]=tensor2lagrangianDescriptors(globalTensor,Tt,countGoodSteps,countNullTensor,lagvelo);
                if Tt==0
                    Lag(isnan(Lag(:)))=1;
                    LagLeft(isnan(LagLeft(:)))=1;%avoid NaN
                    %For some we dont want this FTLE should be zero
                end
                
                %Descriptor Right from tinilag->tfinal
                DescriptorsLs(count,1:2)=[cellid{tini_mat}(:,cell_ini) track_id];
                %We substitute the cellnum by the track_id
                %cellnum{tini_mat}(:,cell_ini)];%we store cell id and cell number
                DescriptorsLs(count,3)=tinilag;
                DescriptorsLs(count,4:6)=coords(:,cell_ini);
                DescriptorsLs(count,7:(6+size(Lag,2)))=Lag;
                t_des_mat=tcurrent+1;
                %Descriptor Left from tfinal->tinilag
                DescriptorsLsLeft(count,1:2)=[cellid{t_des_mat}(:,cellprev) track_id];
                %We substitute the cellnum by the track_id
                %cellnum{t_des_mat}(:,cellprev)];%we store cell id and cell number
                DescriptorsLsLeft(count,3)=tcurrent;
                DescriptorsLsLeft(count,4:6)=XYZ{t_des_mat}(:,cellprev);
                DescriptorsLsLeft(count, 7:(6+size(LagLeft,2)))=LagLeft;
                
                count=count+1;
                
            end%Finish loop trajectories
            
            integrated=DescriptorsLs(:,stepcheck);
            maxT=max(integrated(:))
            minT=min(integrated(:))
            
            if debug_this
                tinilag
                tini_mat
                tfinal
                tcurrent
                Tt
                maxT
                minT
                pause
            end
            %Discarding things that should not be stored
            DescriptorsLs=DescriptorsLs(DescriptorsLs(:,stepcheck)==Tt,:);
            DescriptorsLsLeft=DescriptorsLsLeft(DescriptorsLsLeft(:,stepcheck)==Tt,:);
            DescriptorsLs=DescriptorsLs(DescriptorsLs(:,3)>=tinilag,:);
            DescriptorsLsLeft=DescriptorsLsLeft(DescriptorsLsLeft(:,3)>=tinilag,:);
            DescriptorsLs=DescriptorsLs(DescriptorsLs(:,1)~=0,:);
            DescriptorsLsLeft=DescriptorsLsLeft(DescriptorsLsLeft(:,1)~=0,:);
            DescriptorsL=vertcat(DescriptorsL, DescriptorsLs);
            DescriptorsLLeft=vertcat(DescriptorsLLeft, DescriptorsLsLeft);
        end
        
        'saving lagrange'
        save(rawLagDescriptors, 'DescriptorsL');%forwards incremental overwrite this
        save(rawLagDescriptorsLeft, 'DescriptorsLLeft');
        
        if save_tracking
            dlmwrite(track_file,track_lag,'delimiter',';','precision', 16);
        end
    end
    
    if export2movit
        FixedWindow =0
        Lag2movit
    end
end


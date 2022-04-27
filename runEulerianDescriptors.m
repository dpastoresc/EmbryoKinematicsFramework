% /* --------------------------------------------------------------------------------------
%  * File:    runEulerianDescriptors.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTANT SPATIAL (EULERIAN) KINEMATICS DESCRIPTORS FROM TRACKING

clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% User's config %%%%%%%%%%%%%%%%%%%%%%%%
dNs=[1]%Dataset id in datasetListBioEmergences.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xtag=''
do_process=0 %0 just exports colormaps if selected
export2movit=1 %can cancel exportation
StartOne=1 %first colormap index
topology=0%compute eulerian topology decriptors
strains=0%compute eulerian incremental strains
velocity=1
lagrange=0
    
for dN=dNs
    %%%% Load parameters of Kinematics descriptors computation
    tinilag=-1
    loadParametersKinematics;
    limit=tfinal;
    
    %%%%%%%%%%%%%%%%%%%%%%% User's config %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% READ ALL THE INTERMEDIATE DATA %%%%%%%%%%%%
    if do_process
        'Read Deformation Tensor Data'
        F=load([folder dataset '_' tensorid '.mat']);
        F = struct2cell(F);
        F=F{1};
        if strains
            %E= load_centData(dataset,strainid);
            E=load([folder dataset '_' strainid '.mat']);
            E = struct2cell(E);
            E=E{1};
        end
        %steps=size(F,2);
        %Load lineage info
        %XYZ=load_centData(dataset,'XYZ');
        XYZ=load([folder dataset '_XYZ.mat']);
        XYZ = struct2cell(XYZ);
        XYZ=XYZ{1};
        %XYZ
        %cellid=load_centData(dataset,'cellid');
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
        Velo = load([folder dataset '_' veloid]);%load_centData(dataset, veloid);
        Velo = struct2cell(Velo);
        Velo=Velo{1};
        MeanVelo=load([folder dataset '_' velomeanid]);%load_centData(dataset, velomeanid);
        MeanVelo = struct2cell(MeanVelo);
        MeanVelo=MeanVelo{1};
        RVelo=load([folder dataset '_' veloRid]);%load_centData(dataset, velomeanid);
        RVelo = struct2cell(RVelo);
        RVelo=RVelo{1};
        
        VeloSl = load([folder dataset '_' veloslid]);%load_centData(dataset, veloid);
        VeloSl = struct2cell(VeloSl);
        VeloSl=VeloSl{1};
        MeanVeloSl = load([folder dataset '_' velomeanslid]);%load_centData(dataset, velomeanid);
        MeanVeloSl = struct2cell(MeanVeloSl);
        MeanVeloSl=MeanVeloSl{1};
        N=load([folder dataset '_' neighid]);%load_centData(dataset, neighid);
        N = struct2cell(N);
        N=N{1};
        
        sizeTotal=0;
        for tx=tini:limit
            t=tx+1;
            sizeTotal=sizeTotal+size(XYZ{t},2);
            %Initialize tables
        end
        if topology
            DescriptorsT=double(zeros(sizeTotal,27));
        end
        if strains
            DescriptorsS=double(zeros(sizeTotal,34));
        end
        if velocity
            DescriptorsV=double(zeros(sizeTotal,17));
        end
        
        count=1;
        for t=tini:limit
            'Kinematics processing step'
            t
            
            tb=t+1;
            coords=XYZ{tb};
            cells=size(coords,2);
            
            %Maybe we should remove the ref, although will ignore it anyway
            %starting cell at 2:cells
            for cell=1:cells
                Tensor=zeros(3,3);
                %pos=coord(:,cell)
                Tensor(:)=F{tb}(:,cell);%3x3 matrix of tensor
                
                %topology, deformations descritpors... look at
                %tensor2topologyDescriptors
                if topology
                    % num2str(cellid{step}(:,cell)) ';' num2str(cellnum{step}(:,cell)) '
                    %[cellid{tb}(:,cell) cellnum{tb}(:,cell)]
                    DescriptorsT(count,1:2)=[cellid{tb}(:,cell) cellnum{tb}(:,cell)];%we store cell id and cell number
                    %[s pp pn qp qn e exp_nv c comp_nv rota rota_nv rat]=tensor2topologyDescriptors(Tensor);
                    %Topo=[s pp pn qp qn e exp_nv c comp_nv rota rota_nv rat];
                    
                    %calculate gain in neighbors
                    %N==[mean(dave) std(dave) mean(anglesrel(:)) std(anglesrel(:)) numOfSamples];
                    n=N{tb}(5,cell);
                    %size(N{tb})
                    %pause
                    Topo=tensor2topologyDescriptors(Tensor, Velo{tb}(:,cell),MeanVelo{tb}(:,cell),n);
                    %Topo(2:3)
                    DescriptorsT(count,3)=t;
                    DescriptorsT(count,4:6)=coords(:,cell);
                    DescriptorsT(count,7:(6+size(Topo,2)))=Topo;
                    %size(DescriptorsT(count,:))
                    %size(Topo)
                    %DescriptorsT(count,8:9)
                end
                
                %strain look at tensor2strainDescritpors
                if strains
                    %Instant strain, could be linear version, but this
                    %is ok.s
                    %Maybe compare with the Strain tensor obtained by
                    %benoit's code, should be the same
                    %Et = 1/2 * ( ( Tensor' * Tensor ) - diag(ones(1,3))  )
                    %if sum(Tensor(:))>0
                    %    Et = 1/2 * ( ( Tensor' * Tensor ) - diag(ones(1,3)) );
                    %else
                    %    Et=zeros([3 3 3]);
                    %end
                    E2=zeros(3,3);
                    E2(:) = E{tb}(:,cell);
                    DescriptorsS(count,1:2)=[cellid{tb}(:,cell) cellnum{tb}(:,cell)];
                    Strain=tensor2strainDescriptors(E2);
                    DescriptorsS(count,3)=t;
                    DescriptorsS(count,4:6)=coords(:,cell);
                    DescriptorsS(count,7:(6+size(Strain,2)))=Strain;
                end
                
                %Velocity and neighbours are always retrieved if the
                %flag is on
                if velocity
                    %calculate gain in neighbors
                    %N==[mean(dave) std(dave) mean(anglesrel(:)) std(anglesrel(:)) numOfSamples];
                    n=N{tb}(5,cell);
                    prevneighs=n;
                    cellDiam=N{tb}(1,cell);
                    convergence=N{tb}(4,cell);
                    if(tb>(tini+1))
                        %size(cellmother{tb})
                        %'just for the cell'
                        li=cellmother{tb}(:,cell);
                        if(li>1)
                            %'previous'
                            %tb-1
                            %'now'
                            %tb
                            %size(N{tb-1})
                            prevneighs=N{tb-1}(5,li);
                        end
                    end
                    %size(VeloSl{tb})
                    %size(MeanVeloSl{tb})
                    
                    DescriptorsV(count,1:2)=[cellid{tb}(:,cell) cellnum{tb}(:,cell)];
                    Mig=vf2migrationDescriptors(Velo{tb}(:,cell),MeanVelo{tb}(:,cell),RVelo{tb}(:,cell),VeloSl{tb}(:,cell),MeanVeloSl{tb}(:,cell),n,prevneighs, cellDiam, convergence);
                    DescriptorsV(count,3)=t;
                    DescriptorsV(count,4:6)=coords(:,cell);
                    %                         size(Mig)
                    %                         Mig
                    %                         pause
                    DescriptorsV(count,7:(6+size(Mig,2)))=Mig;
                    %DescriptorsV(count,12:14)
                end
                
                count=count+1;
            end
        end
        
        %save data
        if topology
            'saving topology'
            save(rawTopoDescriptors, 'DescriptorsT');
        end
        if strains
            'saving strain'
            save(rawStrainDescriptors, 'DescriptorsS');
        end
        if velocity
            'saving velocity'
            save(rawVelDescriptors, 'DescriptorsV');
        end
    end
    %export data to BioEmergences Platform
    %See libMovit/Kinematics2movit for details
    %Configure colormap data in libMovit/loadParametersKinMovit
    if export2movit
        Kinematics2movit
    end
end


